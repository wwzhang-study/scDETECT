#' @title Calculation of projection matrix.
#' @description It calculates the projection matrix of data linear model based on
#' design matrix.
#' @param Y_raw A list of matrix of pseudo-bulk data for each cell type.
#' @param Design_matrix the design matrix of all the covariates of the data
#'   model for gene expression.
#' @param z.states.k DE states in each cell type(0: normal; 1: disease).
#' @param factor.to.test A phenotype name, e.g. "disease", or a vector of
#'   contrast terms, e.g. c("disease", "case", "control").
#' @return the projection matrix

#'
M.cal <- function(Y_raw, Design_matrix, z.states.k, factor.to.test){
  M.reduce = list()
  dz.combine.num = nrow(z.states.k)

  for(dz.ix in 1:dz.combine.num){
    if(z.states.k[dz.ix,]==1){
      design.reduce <- Design_matrix
    }else{
      name.index <- grep(factor.to.test, colnames(Design_matrix))
      design.reduce <- Design_matrix[,-name.index]
    }

    ### projection matrix for column spaced defined on design.reduce
    M.reduce[[dz.ix]] <- design.reduce %*% solve(t(design.reduce) %*% design.reduce) %*% t(design.reduce)
  }
  return(M.reduce)
}
