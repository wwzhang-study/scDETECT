#' @title Calculation of marginal probability.
#' @description It calculates marginal probability by summing over all possible conditions.
#' @param a one single probability.
#' @return the marginal probability.

#'
Rsumlog.matrix <- function(a){
  s <- a[,1]
  for( i in 2:dim(a)[2] ){
    s <- Raddlog.matrix(s, a[,i])
  }
  return(s)
}
