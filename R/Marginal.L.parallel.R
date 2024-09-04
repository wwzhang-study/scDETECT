#' @title Calculation of data likelihood
#' @description It calculates observed data likelihood based on the data model and
#' parameter estimation in OLS.
#' @param X the projection matrix of data linear model.
#' @param Y_raw A list of matrix of pseudo-bulk data for each cell type.
#' @param geneNames names of all genes.
#' @param all.z.states all DZ combinations based on independent cell types.
#' @param numCores number of cores for parallel running.
#' @return data likelihood
#' @importFrom Seurat LogNormalize
#' @importFrom stats dnorm
#' @importFrom doParallel registerDoParallel
#' @import parallel

#'
Marginal.L.parallel <- function(X, Y_raw, geneNames, all.z.states, numCores){
  N <- nrow(all.z.states)
  #sample.size <- ncol(Y_raw[[1]])

  foo <- function(i, Y_raw, X, geneNames) {

    gene.num = length(geneNames)
    cell.num = length(Y_raw)
    all.z.states <- expand.grid(rep(list(0:1), cell.num))
    p.ycz = matrix(NA, ncol = cell.num, nrow = gene.num)

    for(cell in 1:cell.num){
      z.gk = all.z.states[i, cell]
      Y = Seurat::LogNormalize(Y_raw[[cell]][geneNames,], scale.factor = 100000)
      Y = as.matrix(Y)
      mu_est = Y %*% t(X[[z.gk + 1]])
      resi.temp = Y - mu_est
      sample.size <- ncol(Y_raw[[cell]])
      sd_est = sqrt(rowMeans((resi.temp - rowMeans(resi.temp))^2) * (1 - 1/(sample.size)))
      p.yczk = dnorm(Y, mu_est, sd_est, log = T)
      p.ycz[, cell] = rowSums(p.yczk)
    }
    rowSums(p.ycz)
  }
  doParallel::registerDoParallel(cores = numCores)
  cl <- parallel::makeCluster(numCores)
  p.ycz.sum <- parallel::parSapply(cl, 1:N, foo, Y_raw, X, geneNames)
  parallel::stopCluster(cl)
  return(p.ycz.sum)
}
