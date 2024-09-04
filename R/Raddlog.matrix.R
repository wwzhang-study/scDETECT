#' @title Getting sum of two probabilities.
#' @description It describes how to calculate sum of two probabilities.
#' @param a one single probability.
#' @param b another single probability.
#' @return sum of two probabilities.


Raddlog.matrix <- function(a,b){
  result <- rep(0, length(a))
  idx1 <- (a>b+200) | is.infinite(b)
  result[idx1] <- a[idx1]

  idx2 <- (b>a+200) | is.infinite(a)
  result[idx2] <- b[idx2]

  idx0 <- !(idx1|idx2)
  result[idx0] <- a[idx0] + log1p(exp(b[idx0]-a[idx0]))
  return(result)
}
