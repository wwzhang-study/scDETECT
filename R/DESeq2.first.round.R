#' @title DESeq2 method for differential analysis
#' @description It conducts statistical test with DESeq2 method
#' for each cell type.
#' @param Y_raw A list of matrix of pseudo-bulk data for each
#'   cell type.
#' @param design covariates representing interested factors to be
#'   tested, such as disease, gender and ethnicity.
#' @return A list including the results from DESeq2 test in each
#'   cell type.
#' @importFrom stats rnbinom
#' @export
#'
#' @examples
#' N <- 40 # simulation a dataset with 40 samples
#' K <- 3 # 3 cell types
#' C <- 10 # 10 cells are simulated in each sample for each cell type
#' P <- 500 # 500 features

#' ### simulate single cell data and transform to pseudo-bulk data
#' design <- data.frame(disease=factor(sample(0:1,size = N,replace=TRUE)))
#' Y = list()
#' for (i in 1:K){
#'   y.pseudo = c()
#'     for (j in 1:N){
#'       y = matrix(rnbinom(P*C, size = 1, mu = 1), nrow = P, byrow = FALSE)
#'       y.pseudo = cbind(y.pseudo, rowSums(y))
#'     }
#'   Y[[i]] = y.pseudo
#'   rownames(Y[[i]]) <- paste0('gene',1:P)
#'}

#' ### conduct statistical test in DESeq2 method
#' DESeq2_res <- DESeq2.first.round(Y,design)
DESeq2.first.round <- function(Y_raw,design){

  ### used to store result
  res <- list()
  ### test for each cell type
  for(cell in 1:length(Y_raw)){
    pseudobulk_Yk = Y_raw[[cell]]

    ### store the DESeq2 result in res[[cell]] list
    res[[cell]] <- getPbyDEseq2(pseudobulk_Yk,design)
  }

  return(res) ### return list res
}
