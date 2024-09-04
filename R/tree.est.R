#' @title Tree estimation
#' @description It estimates the cell type hierarchical tree structure based
#' on test statistics from DESeq2 method.
#' @param tstat matrix of test statistics from DESeq2.
#' @param de.res prior DE state for tree estimation.
#' @param similarityFun Custom function used to calculate similarity between
#'   cell types that used for tree structure estimation.
#' @param tree.type Tree type for inference, default is c('single','full').
#' @return A list containing two tree structures, the one is the simplest
#'   tree structure with only one layer and the other is the tree structure
#'   with multiple layers.
#' @importFrom stats as.dist
#' @importFrom stats cor
#' @importFrom stats hclust
#' @importFrom stats cutree

#'
tree.est <- function(tstat, de.res, similarityFun = NULL, tree.type){

  tree.input <- list()
  cell.num <- ncol(tstat)

  if('single' %in% tree.type){
    tree.input[['single']] <- rbind(rep(1,cell.num), seq(1,cell.num,1))
    colnames(tree.input[['single']]) <- colnames(tstat)
  }

  if('full' %in% tree.type){
    tstat0 <- tstat[rowSums(de.res) > 0, ]
    if(is.null(similarityFun)){
      dist.corr <- as.dist( (1 - cor(tstat0, method = 'pearson' ))/2 )
    }else{
      dist.corr <- as.dist(similarityFun(tstat0))
    }

    hc.tree <- hclust( dist.corr )

    cut.height <- sort(c(0,hc.tree$height),decreasing = T )
    cut.thres <- cut.height

    tree.input[['full']] <- t(cutree(hc.tree, h=cut.thres))
    colnames(tree.input[['full']]) <- colnames(tstat)
  }

  tree.input
}
