#' @title Calculation of prior probs.
#' @description It estimates prior probs of DE for each node of tree structure, and
#' calculates joint prior probs for each DZ combination.
#' @param var.index the tree structure with all nodes numbered sequentially.
#' @param dz.combine all DZ combinations based on estimated tree.
#' @param DESeq2_res the results from DESeq2 test in each cell type.
#' @param fdr cut off of fdr used to define DE state to estimate prior prob.
#' @param p.value cut off of p-value used to define DE state to estimate prior prob.
#' @param p.matrix Prior probability on each node of the tree structure.
#' @return A list containing estimates of prior probs (weights) for each DZ combination
#'   and conditional probability for each node.

#'
W.cal <- function(var.index, dz.combine, DESeq2_res, fdr = 0.05, p.value = NULL, p.matrix = NULL){
  dz.combine <- as.matrix(dz.combine)

  if(is.list(DESeq2_res)){
    cell.num <- length(DESeq2_res)
    gene.num <- nrow(DESeq2_res[[1]])
  }else if(is.matrix(DESeq2_res)){
    cell.num <- ncol(DESeq2_res)
    gene.num <- nrow(DESeq2_res)
  }

  layer.num <- nrow(var.index)
  dz.combine.num <- nrow(dz.combine)
  fdr.all.cell <- matrix(NA,ncol = cell.num,nrow = gene.num )

  if(is.list(DESeq2_res)){
    for(i in 1:cell.num){
      if( is.null(p.value)){
        fdr.all.cell[,i] = DESeq2_res[[i]]$padj
      }else if( is.null(fdr)){
        fdr.all.cell[,i] = DESeq2_res[[i]]$pvalue
      }
    }
  }else if(is.matrix(DESeq2_res)){
    fdr.all.cell = DESeq2_res
  }

  if( is.null(p.value)){
    de.all.cell <- (fdr.all.cell < fdr)*1
  }else if( is.null(fdr)){
    de.all.cell <- (fdr.all.cell < p.value)*1
  }

  if(is.null(p.matrix)){
    p.cond.info <- matrix(NA, ncol = cell.num, nrow = layer.num)
    ### Calculate conditional probability
    for(l in (layer.num -1):1){
      nodes <- unique(var.index[l,])
      node.num <- length(nodes)

      for( node in nodes){

        cell.ix <- which(var.index[l,] == node)
        child.var <- unique(var.index[l+1,cell.ix])

        if(length(child.var) ==1){
          p.cond.info[(l+1),cell.ix] <- 1
        }else{
          de.all.ix <- which(rowSums(de.all.cell[,cell.ix]) > 0)

          for(child in child.var){
            ix <- which(var.index[l+1, ] == child)
            if(length(ix) == 1){
              de.num.tmp <- length(which(de.all.cell[,ix] >0 ) )
              if(de.num.tmp ==0){
                de.num.tmp <- 1
              }
              p.cond.info[l+1, ix ] <- de.num.tmp/ length(de.all.ix)
            }else{
              de.num.tmp <- length(which(rowSums(de.all.cell[,ix]) > 0))
              if(de.num.tmp ==0){
                de.num.tmp <- 1
              }
              p.cond.info[l+1, ix ] <- de.num.tmp/ length(de.all.ix)
            }
          }
        }

        if(l == 1){
          if(length(cell.ix) ==1){
            p.cond.info[l, cell.ix] <- mean(de.all.cell[,cell.ix]>0, na.rm = T)
          }else{
            p.cond.info[l, cell.ix] <- mean(rowSums(de.all.cell[,cell.ix])>0, na.rm = T)
          }
        }

      }
    }
  }else{
    p.cond.info <- p.matrix
  }

  ###### calculate prior probs (weights) for each DZ combination
  prod.keep <- matrix(0, ncol= cell.num, nrow= layer.num )
  for(i in 1:layer.num){
    nodes <- unique(var.index[i,])
    for(node in nodes){
      prod.keep[i, which(var.index[i, ] == node )[1]] <- 1
    }

  }

  weights.all <- rep(NA, dz.combine.num )
  for(dz.combine.ix in 1:dz.combine.num){
    dz.state <- dz.combine[dz.combine.ix, ]
    #cat(dz.state,'\n')
    dz.state.matrix <- matrix(dz.state[c(var.index)], nrow = layer.num, ncol = cell.num, byrow=F)

    dz.state.cond.matrix.1 <- rbind(dz.state.matrix[1,], dz.state.matrix[-layer.num,] )
    dz.state.cond.matrix.2 <- rbind(1 - dz.state.matrix[1,], dz.state.matrix[-layer.num,])

    dz.state.power.1 <- dz.state.matrix * dz.state.cond.matrix.1
    dz.state.power.2 <- (1 - dz.state.matrix) * dz.state.cond.matrix.2

    p.matrix <- (p.cond.info ^ dz.state.power.1 ) * ( (1 - p.cond.info) ^ dz.state.power.2)

    weights.all[dz.combine.ix] <- prod( p.matrix ^ prod.keep)
  }


  res <- list('weight' = weights.all, 'est_prob'= p.cond.info)
  return(res)
}
