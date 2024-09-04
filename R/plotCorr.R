#' @title Showing DE/DM state correlation between cell types
#' @description It generates test statistics for each pair of
#' cell types and calculate corresponding Pearson correlation
#' and odds ratio.
#' @param tstat matrix of test statistics from DESeq2.
#' @param de.state matrix of DE/DM states (1: DE/DM, 0:non-DE/DM)
#' @param tstat.thres threshold of test statistics to define DE/DM,
#' required if de.state not provided.
#' @param p.size point size for scatter plot.
#' @param p.color point color for scatter plot.
#' @param fig.margin figure margin.
#' @param fig.margin.unit unit of figure margin.
#' @param line.type line type in scatter plot.
#' @param line.color line color in scater plot.
#' @return A figure contains scatter plots, Pearson correlation and
#'   odds ration of test statistics for each pair of cell types.
#' @import ggplot2 GGally
#' @importFrom grDevices adjustcolor
#' @importFrom stats cor.test
#' @importFrom stats fisher.test
#' @export

#'
#' @examples
#' tstat.1 <- runif(1000,0,5)
#' tstat.2 <- tstat.1 + rnorm(1000,0,1)
#' tstat.2[tstat.2 < 0] =0
#' tstat.3 <- runif(1000,0,5)

#' tstat.input <- data.frame('cell.1'=tstat.1,
#'                           'cell.2'=tstat.2,
#'                           'cell.3'=tstat.3)

#' plotCorr(tstat = tstat.input, tstat.thres = 2.58)
plotCorr <- function(tstat,
                     de.state=NULL,
                     tstat.thres=NULL,
                     p.size = 0.2,
                     p.color = grDevices::adjustcolor( "black", alpha.f = 0.2),
                     fig.margin = c(1,1,1,1),
                     fig.margin.unit = 'in',
                     line.type = 'dashed',
                     line.color = 'blue'
){
  cell.names <- colnames(tstat)
  if(!is.data.frame(tstat)){
    if(is.matrix(tstat)){
      tstat <- data.frame(tstat)
    }else{
      stop('tstat should be a data frame or matrix')
    }
  }

  dim.tstat <- dim(tstat)

  if(is.null(cell.names)){
    message('No cell names detected, new names created for each cell type. \n')
    cell.names <- paste0('cell.',seq(1,dim(tstat)[2],1))
  }

  if(is.null(de.state) & is.null(tstat.thres)){
    tstat.thres = 2.58
    message('No threshold as input to define DE/DMC. Use tstat = 2.58 as threshold.\n')
  }

  if(!is.null(de.state)){
    message('Detect input of de.state.')
    if(!is.data.frame(de.state)){
      if(is.matrix(de.state)){
        de.state <- data.frame(de.state)
      }else{
        stop('de.state should be a data frame or matrix \n')
      }
    }

    dim.de.state <- dim(de.state)
    if( ! all(dim.tstat == dim.de.state) ){
      stop('dimensions of two inputs do not match\n')
    }

  }else if(is.null(de.state) & !is.null(tstat.thres) ){
    de.state <- matrix(NA, ncol= ncol(tstat), nrow = nrow(tstat))
    if(length(tstat.thres)==1){
      de.state <- (abs(tstat) > tstat.thres)*1
    }else if(length(tstat.thres)==ncol(tstat)){
      for( i in 1:ncol(tstat)){
        de.state[,i] <- (abs(tstat[,i]) > tstat.thres[i])*1
      }
    }else{
      stop('tstat.thres length is not correct: assign one value for all cell types
           or assign different values to each cell type.')
    }

    message('Detect input of tstat.thres; Use tstat.thres to calculate odds ratio.\n')
  }

  colnames(tstat) <- cell.names
  colnames(de.state) <- cell.names
  de.state <- data.frame(de.state)

  cell.num <- dim.tstat[2]

  cell.thres <- rep(NA, cell.num)
  for( cell.ix in 1:cell.num ){
    cell.thres[cell.ix] <- min(abs(tstat)[de.state[,cell.ix]==1,cell.ix])
  }
  message(paste0(c("test statistics threshold for each cell type:\n", round(cell.thres,digits=3) ), collapse = " " ))

  fig.tmp <- GGally::ggpairs(data = tstat,
                             upper = list(continuous = GGally::wrap(my_custom_cor,
                                                            de.info = de.state)),
                             lower = list(continuous = GGally::wrap("points", size = p.size,
                                                            colour = p.color)),
                             axisLabels = 'internal')

  for( i in 2:cell.num){
    for(j in 1:(i-1)){
      fig.tmp[i,j] <- fig.tmp[i,j] +
        ggplot2::geom_vline(xintercept = cell.thres[j], linetype = line.type,
                   colour= line.color ) +
        ggplot2::geom_hline(yintercept = cell.thres[i], linetype = line.type,
                   colour= line.color)
    }
  }

  fig.tmp + ggplot2::theme(plot.margin = ggplot2::unit(fig.margin, fig.margin.unit))

}




my_custom_cor <- function(data, mapping, de.info, color = 'black',
                          sizeRange = c(1,5), ...){

  # get the x and y data to use the other code
  x <-GGally::eval_data_col(data, mapping$x)
  y <-GGally::eval_data_col(data, mapping$y)

  ct <- cor.test(x,y) ### test -log10(pval) correlation

  ct.symbol <- c("***", "**", "*", ".", " ")[(ct$p.value  < c(0.001, 0.01, 0.05, 0.1, 1) )][1]
  r <- unname(ct$estimate)
  rt <- format(r, digits=2)[1]

  ## odds ratio
  de.x <- GGally::eval_data_col(de.info, mapping$x)
  de.y <- GGally::eval_data_col(de.info, mapping$y)

  m.tmp <- cbind(de.x, de.y)

  a <- as.numeric( sum( rowSums(m.tmp) == 2 ) ) + 1
  b <- as.numeric( sum( (de.x - de.y) == 1 )) + 1
  c <- as.numeric( sum( (de.x - de.y) == -1 ) ) + 1
  d <- as.numeric( sum( rowSums(m.tmp) == 0 ) ) + 1

  or.est <- format( d/b/c*a, digits=2)
  or.pval <- fisher.test(x=matrix(c(a,c,b,d),2,2))$p.value

  or.symbol <- c("***", "**", "*", ".", " ")[( or.pval <=  c(0.001, 0.01, 0.05, 0.1, 1) )][1]

  cex <- max(sizeRange)

  GGally::ggally_text(
    label = paste0( 'Corr: ',as.character(rt), ct.symbol,'\n',
                    'OR: ', as.character(or.est), or.symbol),
    mapping = ggplot2::aes(),
    xP = 0.5, yP = 0.5,
    size = 5,
    color = color,
    ...
  ) +
    ggplot2::geom_text(
      ggplot2::aes_string(
        x = 0.8,
        y = 0.8
      ),
      label = " ",
      size = 1,
      color = color,
      ...
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(
        color = color,
        linetype = "longdash"
      ),
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank()
    )



}

