#' @title All DZ combinations of all nodes of the hierarchical tree.
#' @description It creates all DZ combinations based on estimated tree structure.
#' @param tree.input the list of two tree structures, one is the simplest
#'   structure with only one layer and the other is the estimated whole tree
#'   structure with multiple layers.
#' @return A list containing all DZ combinations based on estimated tree, as well
#'   as the tree structure with all nodes numbered sequentially.
#' @importFrom tidyr expand_grid
#' @export
#'
DZ_combination_gen <- function(tree.input ){

  tree.level <- dim(tree.input)[1]
  cell.num   <- dim(tree.input)[2]

  var.index <- matrix(NA, ncol = cell.num, nrow = tree.level)
  var.num <- rep(NA, tree.level)
  for( i in 1:tree.level){
    var.num[i] <- length( unique(tree.input[i,]) )
    var.index[i,] <- tree.input[i,]  + sum(var.num[1:(i-1)])*( as.integer(i>1))
  }

  combine.current.tmp <- NULL

  for(tl.ix in 1:tree.level){
    tl.element.num <- length(unique(var.index[tl.ix,]))
    combine.single.tmp <-  expand.grid( rep(list(0:1), tl.element.num))
    combine.current.tmp <-  as.matrix( tidyr::expand_grid( x1=combine.current.tmp, x2 = combine.single.tmp ) )
    ## filter criteria:
    ## 1. splitting node: 0 -> 0, 1 -> 1 or 0;
    ## 2. single node keep previous result
    ## start from second layer of tree to the end of the tree
    if(tl.ix > 1){
      for( cell in 1:cell.num){
        if(tl.ix == 2){
          ix.rm <-  which(  apply(as.matrix(combine.current.tmp[,var.index[1:tl.ix,cell]]),1,diff) == 1 )
        }else{
          ix.rm <-  which( apply( apply(as.matrix(combine.current.tmp[,var.index[1:tl.ix,cell]]),1,diff), 2, max) == 1 )
        }
        if(length(ix.rm) > 0 ){
          combine.current.tmp <- combine.current.tmp[-ix.rm,]
        }

      }

      # criteria 2
      ## 2. single node keep previous result
      ## strategy: a. check whether up node split; b. if not split, we only keep combination that keep same between the two nodes in two layers
      tl.up.ix <- tl.ix - 1
      nodes.up <- var.index[tl.up.ix, ]
      nodes.current <- var.index[tl.ix,]
      nodes.up.unique <- unique(nodes.up)
      for(nodes.up.ix in nodes.up.unique){
        if( length(unique(nodes.current[nodes.up == nodes.up.ix])) == 1 ){
          ix.rm <- which( apply(combine.current.tmp[ ,  c(nodes.up.ix, unique(nodes.current[nodes.up == nodes.up.ix]) ) ], 1, diff) !=0 )
          if(length(ix.rm) > 0 ){
            combine.current.tmp <- combine.current.tmp[-ix.rm,]
          }
        }


      }
    }


  }

  return(list(node.index = var.index, dz.combine = combine.current.tmp))
}
