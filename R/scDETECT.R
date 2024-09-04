#' @title The main function of scDETECT. It performs scRNA-seq DE analysis
#'    accounting for the cell type correlations.
#' @description It provides posterior probability of whether a feature
#'    is DE in certain cell type given pseudo-bulk data.
#' @param Y_raw A list of matrix of pseudo-bulk data for each cell type,
#'   with rows representing features and columns representing samples.
#' @param design.1 covariates representing interested factors to be tested,
#'   such as disease, gender and ethnicity.
#' @param design.2 other covariates that may also affect expression level
#'   beside the factor of interest.
#' @param factor.to.test A phenotype name, e.g. "disease", or a vector of
#'   contrast terms, e.g. c("disease", "case", "control").
#' @param pval matrix of p-values from DESeq2.
#' @param p.adj matrix of adjusted p-values from DESeq2.
#' @param tree hierarchical tree structure used to account cell type correlation.
#' @param p.matrix.input Prior probability on each node of the tree structure.
#' @param de.state DE/DM state of each feature in each cell type.
#' @param cutoff.tree Cut off used to define DE state to estimate tree, it could
#'   be 'fdr', 'pval' or 'tstat'.
#' @param cutoff.prior.prob Cut off used to define DE state to estimate prior probs
#'   of nodes on tree, it could be 'fdr' or 'pval'.
#' @param similarity.function Custom function used to calculate similarity between
#'   cell types that used for tree structure estimation.
#' @param parallel.core Number of cores for parallel running.
#' @param corr.fig A boolean value, whether to plot corrrelation between cell types
#'   using function plotCorr().
#' @param run.time A boolean value, whether to report running time in seconds.
#' @param tree.type Tree type for inference, default is c('single','full').
#' @return A list containing lists for every cell type. Each list contains results
#'   from DESeq2 method, estimates of conditional probability for each node, prior
#'   probs (weights) for each DZ combination, all DZ combinations based on estimated
#'   tree structure, the hierarchical tree structure and the posterior probability
#'   of each gene in each cell type.
#' @importFrom stats p.adjust
#' @importFrom stats na.omit
#' @export
#' @examples
#' N <- 40 # simulation a dataset with 40 samples
#' K <- 3 # 3 cell types
#' C <- 10 # 10 cells are simulated in each sample for each cell type
#' P <- 500 # 500 features

#' ### simulate single cell data and transform to pseudo-bulk data
#' design <- data.frame(disease=factor(sample(0:1,size = N,replace=TRUE)))
#' Y = list()
#' for (i in 1:K){
#'   y.pseudo = c()
#'     for (j in 1:N){
#'       y = matrix(rnbinom(P*C, size = 1, mu = 1), nrow = P, byrow = FALSE)
#'       y.pseudo = cbind(y.pseudo, rowSums(y))
#'     }
#'   Y[[i]] = y.pseudo
#'   rownames(Y[[i]]) <- paste0('gene',1:P)
#'}

### run scDETECT
#' res <- scDETECT(Y_raw = Y,
#'                 design.1 = design,
#'                 design.2 = NULL,
#'                 factor.to.test = 'disease',
#'                 cutoff.tree = c('tstat',2.58),
#'                 corr.fig = TRUE,
#'                 cutoff.prior.prob = c('pval',0.1) )
scDETECT <- function(Y_raw, # A list composed of pseudo-bulk data for each cell type
                     ## row: genes, col: samples
                     design.1, # covariates representing interested factors to be tested
                     # row: sample, col: covariate
                     design.2 = NULL, # other covariates that may also affect expression level
                     # row: sample, col: covariate
                     factor.to.test = NULL, # covariate to be tested
                     # form1: only name of covariate
                     # form2: covariate name + contrast levels
                     #        (ref level at last)
                     pval = NULL, # independent inference result "p-value" from DESeq2
                     # row: feature, col: cell type
                     p.adj = NULL, # independent inference result "fdr" from DESeq2
                     # row: feature, col: cell type
                     tree = NULL, # tree used to account cell type correlation
                     # col: cell type, row: tree layer
                     # in same row, different numbers represent different nodes
                     # cell types with same number in same row
                     # means they have same internal node at this level.
                     # example:
                     #   1,1,1,1
                     #   1,1,2,2
                     #   1,2,3,4
                     # If assign cell type name for each column,
                     # the name must be the same as output of 'pval'
                     p.matrix.input = NULL, # prior probability on each node of the
                     # tree structure. Only work when tree structure has been specified.
                     # The dimension must same as tree input.
                     de.state = NULL, # de.state of each feature in each cell type
                     # 0: non-DE; 1: DE
                     cutoff.tree = c('fdr', 0.01), # cut off used to define DE
                     # state to estimate tree
                     # could be 'fdr', 'pval' or 'tstat'
                     # default it 'fdr'=0.01
                     cutoff.prior.prob = c('pval', 0.01), # cut off used to
                     # define DE state to
                     # estimate prior prob
                     # could be 'fdr' or 'pval'
                     # default is 'pval'= 0.01
                     similarity.function = NULL, # function used to calculate similarity
                     # between cell types. The input is matrix of test statistic
                     # dimension is : selected gene number * cell number
                     parallel.core = NULL, # integer to specificy number of
                     # cores to use
                     corr.fig = FALSE,     # whether to plot tstat correlation,
                     run.time = TRUE, # whether report running time
                     tree.type = c('single','full') # two tree structures as input
){

  ### Basic info
  gene.num <- dim(Y_raw[[1]])[1]
  cell.num <- length(Y_raw)
  celltypes = names(Y_raw)

  #########################
  ## Step1: run DESeq2 to get test statistics for cell type hierarchy
  ##        and prior probability
  ########################

  ## make design matrix for posterior probability calculation in later steps
  Design_out <- makeDesign_tree(design.1 = design.1, design.2 = NULL,
                                factor.to.test = factor.to.test)
  Design_matrix <- Design_out$design_matrix

  DESeq2_res <- NULL
  time_res <- list()
  if( is.null(pval) & is.null(p.adj) ){
    message('No prior inference information, run DESeq2 for first round inference \n')

    time.tmp1 <- Sys.time()
    DESeq2_res <- DESeq2.first.round(Y_raw = Y_raw,design = design.1)
    time_res[['DESeq']] <- as.numeric(Sys.time()) - as.numeric(time.tmp1)

    ### extract pval,fdr and test statistic information from first round DESeq2 analysis
    pval = p.adj = tstat <- matrix(NA, nrow = gene.num, ncol = cell.num)
    colnames(pval) = colnames(p.adj) = colnames(tstat) <- celltypes
    rownames(pval) = rownames(p.adj) = rownames(tstat) <- rownames(Y_raw[[1]])

    for(i in 1:cell.num ){
      pval[,i] <- DESeq2_res[[i]]$pvalue
      p.adj[,i] <- p.adjust(pval[,i],method = "BH")
      tstat[,i] <- DESeq2_res[[i]]$stat
    }

  }else if( !is.null(pval) & is.null(p.adj) ){
    cell.types <- colnames(pval)
    p.adj <- matrix(NA, nrow = gene.num, ncol = cell.num )
    colnames(p.adj) <- cell.types
    rownames(p.adj) <- rownames(pval)
    ### calculate fdr for each cell type p-value.
    for( i in 1:cell.num ){
      p.adj[,i] <- p.adjust(pval[,i],method='fdr')
    }

  }else if( is.null(pval) & !is.null(p.adj) ){
    stop( 'p-value is needed')
  }

  pval = na.omit(pval)
  p.adj = na.omit(p.adj)
  tstat = na.omit(tstat)

  ### determine DE state for tree estimation if not provided
  de.res <- matrix(NA,ncol = cell.num, nrow = nrow(tstat))

  if(is.null(de.state)){
    # if de.state not provided then
    # use fdr as cutoff
    if(cutoff.tree[1] == 'fdr'){
      if(length(cutoff.tree) == (cell.num + 1) ){
        # in this case, each cell type has different cutoffs
        for(cell.ix in 1:cell.num){
          de.res[,cell.ix] <- (p.adj[,cell.ix] <
                                 as.numeric(cutoff.tree[,(cell.ix+1)]))*1
        }

        if( sum( rowSums(de.res) > 0) < 2 & 'full' %in% tree.type){

          cat('Not enough DE called for tree estimation. \n
              Please use less restrictive cutoff for tree estimation.')
          if(!is.null(DESeq2_res)){
            return(list('DESeq2_res' = DESeq2_res))
          }
          stop('scDETECT stops')
        }

      }else if(length(cutoff.tree) == 2){
        # in this case, all cell types have same cutoffs
        de.res <- (p.adj < as.numeric(cutoff.tree[2])*1)

        if( sum( rowSums(de.res) > 0) < 2 & 'full' %in% tree.type){

          cat('Not enough DE called for tree estimation. \n
              Please use less restrictive cutoff for tree estimation.')
          if(!is.null(DESeq2_res)){
            return(list('DESeq2_res' = DESeq2_res))
          }
          stop('scDETECT stops')
        }

      }else{
        stop("length of cutoff.tree is incorrect.")
      }

    }else if(cutoff.tree[1]=='pval'){
      # use pvalue as cutoff
      if(length(cutoff.tree) == (cell.num + 1) ){
        for(cell.ix in 1:cell.num){
          de.res[,cell.ix] <- (pval[,cell.ix] <
                                 as.numeric(cutoff.tree[,(cell.ix+1)]))*1
        }
      }else if(length(cutoff.tree) == 2){
        de.res <- (pval < as.numeric(cutoff.tree[2]))*1
      }else{
        stop("length of cutoff.tree is incorrect.")
      }
    }else if(cutoff.tree[1]=='tstat'){
      if(length(cutoff.tree) == (cell.num + 1) ){
        for(cell.ix in 1:cell.num){
          de.res[,cell.ix] <- (abs(tstat[,cell.ix]) >
                                 as.numeric(cutoff.tree[,(cell.ix+1)]))*1
        }
      }else if(length(cutoff.tree) == 2){
        de.res <- (abs(tstat) > as.numeric(cutoff.tree[2]))*1
      }else{
        stop("length of cutoff.tree is incorrect.")
      }
    }else{
      stop('Invalid input of cutoff.tree variable')
    }
  }else{
    de.res <- de.state
  }

  ### Step 1.5: plotCorr show correlation between cell types based on first round
  ###           independent test result
  fig.res <- NULL
  if(corr.fig == TRUE){
    message('Generating scatter plot of tstat among cell types \n')
    fig.res <- plotCorr(tstat = data.frame(tstat),
                        de.state = data.frame(de.res))
  }

  ### Step 2: Estimate a tree structure based on tstat
  ###         Features used should be filtered by users
  ###         This includes two modes: 1. user gives features to build up trees
  ###                                  2. Estimate a tree by FDR/pval cutoff

  if( is.null(tree) ){
    tree.input <- tree.est(tstat, de.res, similarityFun = similarity.function, tree.type)
  }else{
    tree.type <- 'custom'
    tree.input <- list()
    tree.input[['custom']] <- tree
  }

  ### Step 3: Estimate prior probabilities based on given tree structure
  ### Step 4: Estimate posterior probability
  ### Input includes:
  ### 1. tree structure
  ### 2. design matrix
  ### 3. observed bulk data
  ### 4. first round pvalue result (could come from other packages)

  if(cutoff.prior.prob[1]=='fdr'){
    cutoff.fdr <- as.numeric(cutoff.prior.prob[2])
    cutoff.pval <- NULL
    inference_input <- p.adj
  }else if(cutoff.prior.prob[1] == 'pval'){
    cutoff.fdr <- NULL
    cutoff.pval <- as.numeric(cutoff.prior.prob[2])
    inference_input <- pval
  }

  tree_res <- list()
  for(tree.ix in tree.type){
    message(paste0('inference with tree: ',tree.ix, '\n'))
    time.tmp2 <- Sys.time()
    tree_res[[tree.ix]]<- post_calc_tree(Design_matrix = Design_matrix,
                                         Y_raw = Y_raw,
                                         tree.input = tree.input[[tree.ix]],
                                         p.matrix.input = p.matrix.input,
                                         DESeq2_res = inference_input,
                                         factor.to.test = factor.to.test,
                                         fdr = cutoff.fdr,
                                         p.value = cutoff.pval,
                                         core.num = parallel.core)
    time_res[[tree.ix]] <- as.numeric(Sys.time()) - as.numeric(time.tmp2)

  }

  if(run.time == TRUE & corr.fig == TRUE){
    all.res <- list('DESeq2_res'=DESeq2_res, 'tree_res'=tree_res, 'fig'=fig.res,
                    'time_used'=time_res)
  }else if(run.time == FALSE & corr.fig == TRUE){
    all.res <- list('DESeq2_res'=DESeq2_res, 'tree_res'=tree_res, 'fig'=fig.res)
  }else if(run.time == TRUE & corr.fig == FALSE){
    all.res <- list('DESeq2_res'=DESeq2_res, 'tree_res'=tree_res, 'time_used'=time_res)
  }else{
    all.res <- list('DESeq2_res'=DESeq2_res, 'tree_res'=tree_res)
  }

  return(all.res)
}
