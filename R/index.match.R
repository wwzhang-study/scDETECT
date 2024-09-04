#' @title Index of DZ combinations
#' @description It matches the index of the DZ combinations based on estimated
#' tree structures among the all DZ combinations based on independent cell types.
#' @param dz.combine DZ combinations based on estimated tree structure.
#' @param all.z.states all DZ combinations based on independent cell types.
#' @return A vector of index.

#'
index.match <- function(dz.combine, all.z.states){
  match.ix <- c()
  cell.num <- ncol(all.z.states)
  col.num <- dim(dz.combine)[2]
  z.cols <- (col.num - cell.num + 1):col.num

  for( i in 1:nrow(dz.combine)){

    for( n in 1:nrow(all.z.states)){
      if(all(dz.combine[i, z.cols] == all.z.states[n,])  ){
        match.ix <- c(match.ix, n)
        break
      }
    }
  }
  return(match.ix)
}
