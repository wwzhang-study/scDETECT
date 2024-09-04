#' @title Posterior probability calculation
#' @description It calculates posterior probability of whether a feature is DE/DM in
#' each cell type based on Bayesian model.
#' @param Design_matrix the design matrix of all the covariates of the data model for
#'   gene expression.
#' @param Y_raw A list of matrix of pseudo-bulk data for each cell type.
#' @param tree.input the list of two tree structures, one is the simplest
#'   structure with only one layer and the other is the estimated whole
#'   tree structure with multiple layers.
#' @param DESeq2_res the results from DESeq2 test in each cell type.
#' @param factor.to.test A phenotype name, e.g. "disease", or a vector of contrast terms,
#'   e.g. c("disease", "case", "control").
#' @param fdr cut off of fdr used to define DE state to estimate prior prob.
#' @param p.value cut off of p-value used to define DE state to estimate prior prob.
#' @param no_prior_info A boolean value, whether there is not prior information.
#' @param p.matrix.input Prior probability on each node of the tree structure.
#' @param core.num number of cores for parallel running.
#' @return A list containing lists for every cell type. Each list contains estimates
#'   of conditional probability for each node, prior probs (weights) for each DZ
#'   combination, all DZ combinations based on estimated tree structure, the hierarchical
#'   tree structure and the posterior probability of each gene in each cell type.
#' @import Seurat doParallel parallel

#'
post_calc_tree <- function(Design_matrix, Y_raw, tree.input, DESeq2_res,
                           factor.to.test, fdr, p.value =NULL, no_prior_info =FALSE,
                           p.matrix.input = NULL,core.num=NULL){

  gene.num <- nrow(DESeq2_res)      ### gene number
  sample.size <- ncol(Y_raw[[1]])       ### sample size
  cell.num <- ncol(tree.input)      ### celltype number

  # create all DZ combination based on estimated tree structure
  dz.combine.tmp <- DZ_combination_gen(tree.input)
  var.index  <- dz.combine.tmp[[1]]
  dz.combine <- dz.combine.tmp[[2]]

  z.states.k <- expand.grid(rep(list(0:1),1))
  z.states <- dz.combine[,(max(var.index)-cell.num+1):(max(var.index))]
  all.z.states <- expand.grid(rep(list(0:1), cell.num))
  colnames(all.z.states) <- colnames(tree.input)

  match.index <- index.match(dz.combine, all.z.states)

  M.reduce <- M.cal(Y_raw, Design_matrix, z.states.k, factor.to.test = 'disease')

  W.all <- W.cal(var.index, dz.combine, DESeq2_res, fdr = fdr,
                 p.value = p.value, p.matrix = p.matrix.input)
  weight.all <- W.all[['weight']]
  est.prob <- W.all[['est_prob']]

  ### different combination numbers between Z's and D
  dz.combine.num <- nrow(dz.combine)
  all.z.num <- nrow(all.z.states)

  ### pp is used to store posterior probability:
  pp <- matrix(NA, ncol = cell.num, nrow = gene.num)

  p.ycz.sum <- matrix(NA, nrow = all.z.num, ncol = gene.num)

  ### log P(Z=1, Y) a vector, in which each element corresponding to a cell type
  p.z1y <- matrix(NA, nrow = gene.num, ncol = cell.num)

  if(is.null(core.num)){
    cores_num <- 1
  }else{
    cores_num <- min(core.num, parallel::detectCores()) #detectCores() - 2
  }

  p.ycz.sum.temp <- Marginal.L.parallel(X = M.reduce, Y_raw = Y_raw,
                                        geneNames = rownames(DESeq2_res),
                                        all.z.states, numCores = cores_num)

  p.ycz.sum <- t(p.ycz.sum.temp)

  if(no_prior_info == TRUE){
    p.yzd <- t(p.ycz.sum)
    p.y <- Rsumlog.matrix(p.yzd)

    for( i in 1:cell.num){
      p.z1y[,i] <- Rsumlog.matrix(p.yzd[,all.z.states[,i]==1])
    }

  }else{

    p.ycz.sum.dz.states <- p.ycz.sum[match.index, ]

    p.yzd <- t(p.ycz.sum.dz.states + log(weight.all+0.1))

    p.y <- Rsumlog.matrix(p.yzd)

    for( i in 1:cell.num){
      p.z1y[,i] <- Rsumlog.matrix(p.yzd[,z.states[,i]==1 ])
    }

  }

  pp <- exp(p.z1y - p.y)
  colnames(pp) <- colnames(tree.input)
  rownames(pp) <- rownames(DESeq2_res)

  res.all <- list('est_prob'= est.prob, 'weight'=weight.all, 'dz_combine'=dz.combine,
                  'tree_structure'=tree.input, 'node_index'=var.index ,'pp'=pp)
  return(res.all)

}
