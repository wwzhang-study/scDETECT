#' @title DE test for each cell type using DESeq2 on pseudo-bulk.
#' @description Apply DESeq2 to call DE between two groups for each cell type.
#' @param pseudobulk_Yk a pseudo-bulk expression matrix for cell type k.
#' @param design covariates representing interested factors to be
#' tested, such as disease, gender and ethnicity.
#' @return A list including the results from DESeq2 test in each cell type.
#' @importFrom DESeq2 DESeqDataSetFromMatrix
#' @importFrom DESeq2 DESeq
#' @importFrom DESeq2 results
#' @importFrom stats rnbinom
#' @importFrom S4Vectors DataFrame
#' @export
#'
#' @examples
#' N <- 40 # simulation a dataset with 40 samples
#' K <- 3 # 3 cell types
#' C <- 10 # 10 cells are simulated in each sample for each cell type
#' P <- 500 # 500 features

#' ### simulate single cell data and transform to pseudo-bulk data
#' design <- data.frame(disease=factor(sample(0:1,size = N,replace=TRUE)))
#' Y_raw = list()
#' for (i in 1:K){
#'   y.pseudo = c()
#'     for (j in 1:N){
#'       y = matrix(rnbinom(P*C, size = 1, mu = 1), nrow = P, byrow = FALSE)
#'       y.pseudo = cbind(y.pseudo, rowSums(y))
#'     }
#'   Y_raw[[i]] = y.pseudo
#'   rownames(Y_raw[[i]]) <- paste0('gene',1:P)
#'}

#' ### conduct statistical test in DESeq2 method
#' res <- list()
#' for(i in 1:length(Y_raw)){
#'   pseudobulk_Yk = Y_raw[[i]]
#'   res[[i]] <- getPbyDEseq2(pseudobulk_Yk,design)
#' }
#'
getPbyDEseq2 <- function(pseudobulk_Yk,design){
  mycond = factor(design[,1])
  dds = DESeqDataSetFromMatrix(pseudobulk_Yk, DataFrame(mycond), ~mycond)
  dds <- DESeq(dds)
  res = results(dds)
  return(res)
}
