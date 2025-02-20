#####################################################################
# Package: Tessa
# Version: 1.0.0
# Date : 2025-2-12
######################################################################

##########################################################
#       Tessa construct object functions and preprocessing related functions   #
##########################################################
#' @title Each Tessa object has a number of slots which store information. Key slots to access
#' are listed below.
#' @slot gene_expression The gene expression matrix with genes in rows
#' and spots in columns.
#' @slot meta_df The location coorinates and pseudotime aframe
#' @slot covariates The covariate design matrix modeling gene expression.
#' #! try covariates in
#' @slot kernel The kernel matrix to capture spatial correlation, trajectory correlation and error term
#' @slot bandwidth The bandwidth parameter of Gaussian kernel.
#' @slot normalized If TRUE, gene_expression is normalized. Otherwise, gene_expression is raw count
#' @slot signature_genes The vector of gene names, which are used to estimate pseudotime. This is for LOO double dipping correction step
#' @slot result The summary of uTSVGs, TVGs and SVGs.
#' @export
setClass("Tessa", slots = list(
  gene_expression = "matrix",
  meta_df = "matrix",
  covariates = "ANY",
  kernel = "list",
  bandwidth = "numeric",
  normalized = "logical",
  signature_genes = "character",
  result = "list"
))

setMethod("print", "Tessa", function(x) {
  cat("An object of class Tessa:\n")
  cat("Number of genes:", nrow(x@gene_expression), "\n")
  cat("Number of cells:", ncol(x@gene_expression), "\n")
})

# setMethod("show", "Tessa", function(x) {
#   cat("An object of class Tessa:\n")
#   cat("Number of genes:", nrow(x@gene_expression), "\n")
#   cat("Number of cells:", ncol(x@gene_expression), "\n")
# })
#



#' @title Create the Tessa object
#' @param counts The gene expression matrix of dimension n x G,
#' with genes in rows and spots in columns.
#' @param normalized If TRUE, gene_expression is normalized. Otherwise, gene_expression is raw count
#' @param meta_df The n x 3 dimensional dataframe contains location coordinates and pseudotime.
#' The rownames of meta_df should match the colnames of counts matrix.
#' The colnames of meta_df should be x, y and t
#' @param covariates The covariate (if any) design matrix of dimension n x p,
#' modeling gene expression.
#' The rownames of covariates matrix should match the colnames of counts matrix.
#' @return returns Tessa object.
#' @export
CreateTessaObject <- function(counts, normalized = FALSE, meta_df, covariates = NULL, signature_genes = NULL){
  if (!is.matrix(counts)) {
    tryCatch({
      counts <- as.matrix(counts)
    }, error = function(e) {
      stop('\'counts\' could not be converted to matrix. Please check that \'counts\' is coercible to matrix, such as a matrix, dgCmatrix, or data.frame.')
    })
  }

  meta_df = as.matrix(meta_df)
  genes <- rownames(counts)
  if(is.null(genes)) {
    stop('\'rownames(counts)\' is null!')
  }

  barcodes <- colnames(counts)
  if(is.null(barcodes)) {
    stop('\'colnames(counts)\' is null!')
  }

  spots <- intersect(barcodes, rownames(meta_df))
  meta_df.use <- meta_df[match(spots, rownames(meta_df)),]
  counts.use <- counts[,match(spots, colnames(counts)),drop =FALSE]

  if(!is.null(covariates)){
    covariates.use <- covariates[,match(spots, rownames(covariates))]
  }else{
    covariates.use <- NULL
  }

  if(is.null(signature_genes)){
    stop('\'signature_genes\' is null!')
  }
  signature_genes.use <- rownames(counts.use)[na.omit(match(signature_genes,rownames(counts.use) ))]

  object <- methods::new(
    Class = "Tessa",
    gene_expression = counts.use,
    meta_df = meta_df.use,
    covariates = covariates.use,
    normalized = normalized,
    signature_genes = signature_genes.use
  )
  return(object)
}


#' Anscombe variance stabilizing transformation: NB
#' @param counts gene expression count matrix
#' @param sv normalization parameter
NormalizeVST <- function(counts, sv = 1) {
  varx = apply(counts, 1, var)
  meanx = apply(counts, 1, mean)
  phi = coef(nls(varx ~ meanx + phi * meanx^2, start = list(phi = sv)))

  ## regress out log total counts
  norm_counts <- log(counts + 1/(2 * phi))
  total_counts <- apply(counts, 2, sum)
  res_norm_counts <- t(apply(norm_counts, 1, function(x){resid(lm(x ~ log(total_counts)))} ))

  return(res_norm_counts)
}# end func


#' @title MinMax Normalization
MinMax_Normalize <- function(mtx){
  min_value <- min(mtx, na.rm = TRUE)
  max_value <- max(mtx, na.rm = TRUE)
  normalized_mtx <- (mtx - min_value) / (max_value - min_value)
  normalized_mtx
}

#' @title Preprocess: normalize count matrix, filter out genes and spots
#' @param object Tessa object
#' @param spot.threshold A number to filter out spots with less than spot.threshold gene counts
#' @param gene.threshold A number to filter out genes expressed less than gene.threshold percent of spots
#' @param normalized  If TRUE, object@gene_expression is normalized, skip normalization step.
#' @import stringr
#' @importFrom stringr str_subset
#' @return Returns Tessa object
#' @export
data_preprocess <- function(object, spot.threshold = 10, gene.threshold = 0.1,
                            normalized = NULL){
  counts <- object@gene_expression
  meta_df <- object@meta_df

  if(!is.null(normalized)){
    object@normalized = normalized
  }
  ## skip if normalization is done
  if (object@normalized) {
    counts.normalized <- counts
    gene.threshold <- -Inf
    spot.threshold <- -Inf
  }

  ## Scale locations and pseudotime
  t_vars = str_subset( colnames(meta_df),'lineage')
  meta_df[,c('x','y')] <- MinMax_Normalize(meta_df[,c('x','y')])
  meta_df[,t_vars ] <- MinMax_Normalize(meta_df[,t_vars ] )
  cat('The location and pseudotime have been scaled.\n')

  ## Spot QC
  spots.use <- which(colSums(counts) >= spot.threshold)
  counts.use <- counts[, spots.use,drop = FALSE]
  meta_df <- meta_df[spots.use,]
  numSpots.removed <- ncol(counts) - ncol(counts.use)
  cat(paste(numSpots.removed, 'spots with gene expression less than', spot.threshold, 'have been removed. \n'))

  ## Gene QC
  # Filter out mt genes
  if (length(grep("mt-", rownames(counts.use), ignore.case = T)) > 0) {
    mt_gene_list <- grep("mt-", rownames(counts.use), ignore.case = T)
    counts.use <- counts.use[-mt_gene_list,,drop= FALSE]
    cat(paste(length(mt_gene_list), 'mt genes have been removed. \n'))
  }
  # Filter out low expressed genes
  # Gene expression rate
  num_gene <- nrow(counts.use)
  ExpRate <- apply(counts.use, MARGIN = 1, FUN = function(data.vector){
    non.zero <- sum(data.vector != 0)
    return(non.zero / length(data.vector))
  })
  gene.use <- intersect(names(ExpRate[ExpRate >= gene.threshold]), rownames(counts.use))
  counts.use <- counts.use[match(gene.use, rownames(counts.use)),]
  numGenes.removed <- num_gene - nrow(counts.use)
  cat(paste(numGenes.removed, 'genes with gene expression rate less than', gene.threshold, 'have been removed. \n'))

  ## Normalize count data
  if (!object@normalized ) {
    cat('Normalizing the count data...\n')
    counts.normalized <- NormalizeVST(counts.use,sv = 1)
    object@normalized <- TRUE
  }

  object@gene_expression <- as.matrix(counts.normalized)
  object@meta_df <- meta_df
  object@signature_genes <- rownames(object@gene_expression)[na.omit(match(object@signature_genes,rownames(object@gene_expression) ))]
  return(object)
}



#' @title Normal the matrix by min max value
#' @param M numeric data matrix
normalize_matrix <- function(M) {
  min_val <- min(M); max_val <- max(M)
  # Apply min-max normalization
  M_normalized <- (M - min_val) / (max_val - min_val)
  return(M_normalized)
}

#' @title Build the Gaussian kernel matrix with normalized distance matrix
#' @param X  N by D numeric data matrix
#' @param sigma Positive scalar that specifies the bandwidth of the Gaussian kernel
gaussian_distance_normalized_kernel = function(X = NULL,sigma = NULL){
  exp(-1 * normalize_matrix(as.matrix((dist(X))^2)) /sigma)
}



#' @title Build the Gaussian kernel matrix for one lineage
#' @description
#' Build the Gaussian kernel matrix to model correlation between spots.
#' The Gaussian kernel bandwidth will be selected automatically:
#' for datasets with 5000 spots or fewer,
#' non-parametric Sheather-Jones' bandwidths are computed for each gene,
#' and the median of these values is used as the common bandwidth;
#' for datasets with more than 5000 spots, Silverman's bandwidths are computed for each gene,
#' and the median of these gene-specific bandwidths is used as the common bandwidth.
#' @param object The Tessa object.
#' @param lineage The name of the lineage to use
#' @return The Tessa object with the Gaussian kernel matrix.
#' @export
build_kernelMatrix_lineage <- function(object,lineage='lineage1'){
  meta_df <- na.omit(object@meta_df[,c('x','y',lineage)])
  counts <- object@gene_expression
  counts <- counts[,colnames(counts)[match(rownames(meta_df), colnames(counts))]]
  n <- nrow(counts)
  m <- ncol(counts)
  # print(n)
  # print(m)

  if (n > 5000) {
    ## Bandwidth selection by Silverman's
    bw_vector <- apply(counts, MARGIN = 1, stats::bw.nrd0)
    bw <- median(na.omit(bw_vector))
    object@bandwidth[lineage] <- bw
  } else {
    ## Bandwidth selection by SJ's
    bw_vector <- apply(counts, MARGIN = 1, stats::bw.SJ)
    bw <- median(na.omit(bw_vector))
    object@bandwidth[lineage] <- bw
  }

  ## Gaussian kernel, with min-max normalized distance
  KK_S <- gaussian_distance_normalized_kernel(X = meta_df[,c('x','y')], sigma = bw) + diag(rep(1e-9, m ))
  KK_T <- gaussian_distance_normalized_kernel(X = meta_df[,lineage], sigma = bw)  + diag(rep(1e-9, m ))
  object@kernel[[lineage]] <- list(kernel_S = KK_S, kernel_T = KK_T, kernel_error = diag(rep(1,ncol(counts))))
  return(object)
}



#' @title Build the Gaussian kernel matrix
#' @param object The Tessa object.
#' @return The Tessa object with the Gaussian kernel matrix.
#' @export
build_kernelMatrix <- function(object){
  t_vars <- str_subset( colnames(meta_df),'lineage')
  for(t_var in t_vars ){
    object <- build_kernelMatrix_lineage(object, lineage = t_var )
  }
  object
}

#' @title Overall Test for uTSVG
#' @param Y The numeric matrix with n spots (row) and k genes (column)
#' @param covariates Covariates X as fixed effect
#' @param kernel_mat_list A list of kernel matrix, as random effect. kernel_S, kernel_T and kernel_error
#' @import pracma CompQuadForm
#' @importFrom pracma pinv
#' @importFrom CompQuadForm davies liu
Test1 = function(Y, covariates = NULL, kernel_mat_list, acc = 1e-7){
  # print('Test1')
  # print(str(Y))
  n = nrow(Y)
  ## setup covariate
  if (is.null(covariates)){
    covariates <- matrix(1, n, 1)
  }
  if(!is.matrix(covariates)){
    covariates <- as.matrix(covariates)
  }
  num_cov <- ncol(covariates)
  ## V = K_s + K_t + K_error
  V = matrix(0,n,n)
  for(i in 1:length(kernel_mat_list)){
    V <- V + kernel_mat_list[[i]]
  }
  # SKAT_RL for overall test
  Xdagger <- pinv(covariates) # X+ = Xdagger = (t(X)%*%X)^(-1)%*%t(X), pseudoinverse of covariates
  PVP <- V - (covariates%*%Xdagger)%*%V
  PVP <- PVP - PVP %*% t(Xdagger) %*% t(covariates)

  ## eigen
  zeros_threshold <- 1e-4
  # phis <- getEigenValues(PVP)
  # print('eigen')
  phis <- Re(eigen(PVP, only.values = TRUE)$values)
  phis[phis<zeros_threshold] <- 0.0
  ## calculate k = rank(PVP)
  k <- sum(phis > 0)
  ## calculate q = dim(ker(PVP) & col(P))
  B <- cbind(V, covariates)
  #!may need a new svd, or deal with B
  q <- n - sum(baseSVD(B) > zeros_threshold)
  # q <- n - sum(svd(B)$d > zeros_threshold)
  if(q==0){q <- 1}
  ## calculate score statistics
  ### calculate nominators
  nominators <- apply((PVP%*%Y) * Y, 2, sum)
  ### calculate dnonimators
  denonimators <- apply((Y - (covariates %*% Xdagger %*% Y)) * Y, 2, sum)
  scores <- nominators / denonimators * (n - num_cov)

  ## define pvlaues
  # print('davies p values')
  pvals <- rep(0,length(scores))
  for(i in seq_len(length(scores))){
    lambda <- c(phis[1:k] - scores[i] / (n - num_cov), ones(1, q) * - scores[i] / (n - num_cov))
    pvals[i] <- davies(0, lambda, acc=acc)$Qq
    ## to aviod exactaly zeros using liu's method
    if(pvals[i]<0){
      pvals[i] <- liu(0, lambda)
    }# end fi
    #if(verbose) {cat(paste("SPARK.SCORE::Davies pvalue 1 = ", pvals[i],"\n"))}
  }
  res_test <- data.frame(geneid = colnames(Y),
                         pvs = pvals,
                         pvs_adj = p.adjust(pvals, method = 'BY') )
  res_test
}


#' @title Run Test1 for all genes
#' @param object The TESSA object
#' @param lineage Name of the lineage
#' @param acc If pvalue has too many zeros, decrease this number to adjust the lower limit of pvalue
#' @return The TESSA object
#' @export
run_Test1_lineage = function(object, lineage = 'lineage1',acc = 1e-7 ){
  # print('run_Test1_lineage')
  # print(lineage)
  meta_df <- na.omit(object@meta_df[,c('x','y',lineage)])
  Y <- object@gene_expression[,colnames(object@gene_expression)[match(rownames(meta_df), colnames(object@gene_expression))]]
  if(!is.null(object@covariates)){
      covariates <- object@covariates[match(rownames(meta_df), colnames(counts))]
  }else{
      covariates <- object@covariates
  }
  # print(str(object@kernel[[lineage]]))
  res_test <- Test1(Y = t(Y), covariates = covariates, kernel_mat_list = object@kernel[[lineage]], acc = acc)
  object@result[[lineage]] <- list(Test1 = res_test )
  object
}


#' @title Run Test1 for one gene, one lineage
#' @param object The TESSA object
#' @param gene The gene to do Test1
#' @param lineage The name of lineage
#' @return pvalue of test1
#' @import slingshot irlba SingleCellExperiment S4Vectors
#' @importFrom irlba prcomp_irlba
#' @importFrom slingshot slingshot
#' @importFrom SingleCellExperiment SingleCellExperiment colData
#' @importFrom S4Vectors SimpleList
#' @export
run_Test1_lineage_LOO_pergene = function(object, gene, lineage = 'lineage1', acc = 1e-7,npcs = 30 ){
    meta_df <- na.omit(object@meta_df[,c('x','y',lineage)])
    Y <- object@gene_expression[,colnames(object@gene_expression)[match(rownames(meta_df), colnames(object@gene_expression))]]

    if(!object@normalized) {
      stop('please normalize the data first')
    }

    ## estimate pseudotime
    signature_genes <- setdiff(object@signature_genes, gene)
    embedding <- prcomp_irlba(t(Y[signature_genes,]), scale. = TRUE, n = npcs)$x
    sim <- SingleCellExperiment(assays = Y,
                                reducedDims = SimpleList(PCA = embedding),
                                colData = DataFrame(clusterlabel = rep("0", ncol(Y)) )
                                )
    sim  <- slingshot(sim, clusterLabels = 'clusterlabel', reducedDim = 'PCA',start.clus="0" )
    loc_df[,'lineage1'] = sim@colData@listData$slingPseudotime_1


  ## replaced by new pseudotime
  object <- CreateTessaObject(counts = Y[gene,,drop =FALSE ], meta_df = loc_df,
                              signature_genes =  signature_genes,
                              covariates = object@covariates, normalized = object@normalized )
  object <- build_kernelMatrix(Tessa.obj)
  object <- run_Test1_lineage(object, lineage = 'lineage1',acc = acc )
  pv <- object@result$lineage1$Test1$pvs[1]
  pv
}

#' @title Run Test1 for all genes, all lineages
#' @param object The TESSA object
#' @param LOO If TRUE, do LOO double dipping correction
#' @param acc If pvalue has too many zeros, decrease this number to adjust the lower limit of pvalue
#' @return The TESSA object
#' @import pracma CompQuadForm
#' @importFrom pracma pinv
#' @importFrom CompQuadForm davies liu
#' @export
run_Test1 = function(object, LOO = FALSE, LOO_pv_threshold = 0.05, npcs = 30, acc = 1e-15 ){ #, kernel_names = c('kernel_S','kernel_T', 'kernel_error')
  t_vars = str_subset( colnames(meta_df),'lineage')
  for(lineage in t_vars){
    object <- run_Test1_lineage(object, lineage = lineage,acc = acc )
    if(LOO){
      test_res <- object@result[[lineage]]$Test1
      DP_genes <- intersect(object@signature_genes, test_res$geneid[test_res$pvs_adj < LOO_pv_threshold])
      #! test
      DP_genes_pvs <- unlist(lapply(DP_genes,function(gene){
         run_Test1_lineage_LOO_pergene(object = object, gene = gene, lineage = lienage, npcs = npcs)
        }))
      test_res$pvs_LOO <- test_res$pvs
      test_res$pvs_LOO[match(DP_genes,test_res$geneid )] <- DP_genes_pvs
      test_res$pvs_adj_LOO <- p.adjust(test_res$pvs_LOO, method = 'BY')
      object@result[[lineage]]$Test1 <- test_res
    }
  }
  object
}

#' @title a faster version of tr(M1 %*% M2)
TT <-
  function(M1, M2)
  {
    nn <- nrow(M1)
    S <- c()
    for (itt in 1 : nn)
    {
      S[itt] <- sum(M1[itt, ]*M2[, itt])
    }
    trace <- sum(S)
    return(trace)
  }


#TODO: try gatson, b/c it is the fastest by
## https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009659
#' @title Individual Test for TVG and SVG
#' @param Y The numeric matrix with n spots (row) and k genes (column)
#' @param covariates Covariates X as fixed effect
#' @param kernel_mat_list A list of kernel matrix, as random effect. kernel_S, kernel_T and kernel_error
#' @import MM4LMM
#' @importFrom MM4LMM MMEst
Test2 <- function(object, gene, K_test, lineage = 'lineage1', cores = 1){
  KList =  object@kernel[[lineage]][setdiff(names(object@kernel[[lineage]]), K_test)]
  K_alt = object@kernel[[lineage]][[K_test]]
  Y = object@gene_expression
  Y = Y[gene,colnames(Y)[match(rownames(K_alt),colnames(Y))]]
  n = length(Y)
  if(!is.null(object@covariates )){
    X <- object@covariates
  }else{
    X <- matrix(1, n, 1)
  }

  ## parameter estimation
  ResultVA <- MMEst(Y=Y , Cofactor = X, VarList = KList, NbCores=cores)
  ## if directly use estimation (may try extract info from algo )
  V_l <-  0
  for(i in seq_along(KList)){
    V_l <- V_l +  ResultVA$NullModel$Sigma2[i] * KList[[i]]
  }
  V_l_inv <- invert(V_l)
  P_l <- V_l_inv - V_l_inv%*%X%*% invert(t(X)%*%V_l_inv%*%X) %*%t(X)%*%V_l_inv

  # beta <- ResultVA$NullModel$Beta
  # Q <- t(Y-X%*%beta)%*%V_l_inv%*%K_alt%*%V_l_inv%*%(Y-X%*%beta)/2
  Q <- (t(Y) %*% P_l %*% K_alt %*% P_l%*% Y )/2
  e <- TT(P_l, K_alt)/2

  num_VC = length(KList) #number of VC in Null model
  Ill_values <- c()
  for(i in seq_along(KList)){
    for(j in seq_along(KList)){
      Ill_values <- c(Ill_values, TT(P_l %*% KList[[i]], P_l %*% KList[[j]]  ))
    }
  }

  I_l_l <- matrix(Ill_values ,nrow = num_VC ,ncol = num_VC, byrow = TRUE )/2
  Il_l <- matrix(unlist(lapply(KList, function(K){TT(P_l %*% K_alt, P_l %*% K )}))
                 ,nrow = 1 ,ncol = num_VC, byrow = TRUE)/2
  Ill <- TT(P_l %*% K_alt, P_l %*% K_alt )/2
  Itt <- Ill-Il_l%*%pinv(I_l_l)%*%t(Il_l) #Ill_tilde
  k <- Itt/e/2; v=2*e^2/Itt
  pvalue <- pchisq(Q/k, df=v, lower.tail=F)

  out <- list(gene = gene,VCs=ResultVA$NullModel$Sigma2, #Beta=beta, #restricted.logLik=lR, fisher.info=Iss/2, itr_num = ResultVA$NullModel$NbIt
                 Score=Q, df=v, scale=k, p.value=pvalue)
  return(out)
}


#' @title Run Test2 for selected genes for one lineage
#' @param genes The genes to test for Test2
run_Test2_lineage = function(object, lineage = 'lineage1', genes = NULL, LOO = FALSE, DP_genes = NULL){
  meta_df <- na.omit(object@meta_df[,c('x','y',lineage)])
  Y <- object@gene_expression[,colnames(object@gene_expression)[match(rownames(meta_df), colnames(object@gene_expression))]]

  res_list <- list()
  for(K_test_name in c( 'kernel_S', 'kernel_T')){
    res <- lapply(genes, function(gene){

      if(LOO & K_test_name == 'kernel_T' & gene %in% DP_genes){
        ## get a new t for DP genes
          signature_genes <- setdiff(object@signature_genes, gene)
          embedding <- prcomp_irlba(t(Y[signature_genes,]), scale. = TRUE, n = npcs)$x
          ##'slingshot'
          sim <- SingleCellExperiment(assays = Y,
                                      reducedDims = SimpleList(PCA = embedding),
                                      colData = DataFrame(clusterlabel = rep("0", ncol(Y)) )
          )
          sim  <- slingshot(sim, clusterLabels = 'clusterlabel', reducedDim = 'PCA',start.clus="0" )
          loc_df[,'lineage1'] = sim@colData@listData$slingPseudotime_1

          object <- CreateTessaObject(counts = Y[gene,,drop =FALSE ], meta_df = loc_df,
                                      signature_genes =  signature_genes,
                                      covariates = object@covariates, normalized = object@normalized )
          object <- build_kernelMatrix(Tessa.obj)
          Test2(object = object, gene = gene, K_test = K_test_name ,lineage = 'lineage1')
      }else{
          Test2(object = object, gene = gene, K_test = K_test_name ,lineage = lineage)
      }
    })
    res_list[[K_test_name]] <- unlist(lapply(res, function(x){x$p.value}))
  }
  res_df <- res_list %>% bind_cols() %>% mutate(geneid = genes) %>%
    rename('kernel_T' = 'TVG_pvs', 'kernel_S' = 'SVG_pvs') %>%
    mutate(TVG_pvs_adj = p.adjust(TVG_pvs, method = 'BY'),
           SVG_pvs_adj = p.adjust(SVG_pvs, method = 'BY'))
  object@result[[lineage]][['Test2']] <- res_df[,c('geneid', 'TVG_pvs', 'TVG_pvs_adj', 'SVG_pvs', 'SVG_pvs_adj')]
  object
}

#' @title Run Test1 for all lineages
#' @param pv_threshold The pvalue threshold for test1, to get the uTSVGs that are significant in Test1
#' @export
run_Test2 = function(object, genes = NULL, pv_threshold = 0.05, LOO = FALSE){

  if(LOO){
    DP_genes <- intersect(genes, object@signature_genes)
  }else{
    DP_genes = NULL
  }

  t_vars = str_subset( colnames(meta_df),'lineage')
  for(lineage in t_vars){
    if(is.null(genes)){
      test1_res <- object@result[[lineage]]$Test1
      if('pvs_adj_LOO' %in% colnames(test1_res)){
        genes <- test1_res$geneid[test1_res$pvs_adj_LOO < pv_threshold ]
      }else{
        genes <- test1_res$geneid[test1_res$pvs_adj < pv_threshold ]
      }
    }
    object <- run_Test2_lineage(object, lineage = lineage, genes = genes, LOO = LOO, DP_genes = DP_genes)
  }
  object

}



