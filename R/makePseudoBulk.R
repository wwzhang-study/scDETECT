#' @title Transforming single cell data to pseudo-bulk data.
#' @description It gets summation of single cell data in each
#' sample for each cell type.
#' @param Yk single cell data from cell type k.
#' @param anno_k information of cells in the kth cell type.
#' @return the pseudo-bulk data
#' @importFrom stats rnbinom
#' @export

#'
#' @examples
#' C <- 1000 # 1000 cells are simulated
#' P <- 500 # 500 features

#' y = matrix(rnbinom(P*C, size = 1, mu = 1), nrow = P, byrow = FALSE)
#' cell_info = data.frame(disease=sample(c("control","case"),size = C,
#'    replace=TRUE), majorType=sample(c("CD4","CD8","Mono"),
#'    size = C,replace=TRUE), sampleID=sample(c("sample1","sample2","sample3",
#'    "sample4","sample5"),size = C,replace=TRUE))
#' ix_k = cell_info$majorType == "CD4"
#' anno_k = cell_info[ix_k,]
#' Yk = as.matrix(y[,ix_k])

#' makePseudoBulk(Yk,anno_k)
#'
makePseudoBulk <- function(Yk,anno_k){
  ID = anno_k$sampleID
  allID = names(table(ID))
  result = matrix(0, nrow = nrow(Yk), ncol = length(allID))
  for(i in 1:length(allID)) {
    ix = ID == allID[i]
    if(sum(ix)!=1){
      result[,i] = rowSums(Yk[,ix])
    }
    else{
      result[,i] = Yk[,ix]
    }
  }
  colnames(result) = allID

  return(result)
}
