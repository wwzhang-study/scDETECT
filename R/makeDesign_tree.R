#' @title Making design matrix
#' @description It creates design matrix of covariates of the data model
#' for gene expression.
#' @param design.1 covariates representing interested factors to be tested,
#'   such as disease, gender and ethnicity.
#' @param design.2 other covariates that may also affect expression level
#'   beside the factor of interest.
#' @param factor.to.test A phenotype name, e.g. "disease", or a vector of
#'   contrast terms, e.g. c("disease", "case", "control").
#' @return A list containing a design matrix, all covariates of the data
#'   model and their names, the formula of the linear model.
#' @importFrom stats relevel
#' @importFrom stats model.matrix
#' @importFrom stats as.formula
#' @export
#'
#' @examples
#' design.1 = data.frame(disease = c(0,1,1,0,0,1,0,0,1,0))
#' design.2 = data.frame(gender = c(1,0,1,0,0,1,1,0,1,1))
#' makeDesign_tree(design.1 = design.1, design.2 = design.2,
#'                              factor.to.test = 'disease')
#'
makeDesign_tree <- function(design.1, design.2, factor.to.test){
  covariate.name.1 <- colnames(design.1)
  covariate.name.2 <- colnames(design.2)

  if(ncol(design.1) > 1){
    if(length(factor.to.test) == 1){
      # In this situation, user wants to test the covariate globally if
      # factor.to.test contains multiple levels.

      # put the covariate to last column for convenience
      covariate.name.1 <- c(covariate.name.1[covariate.name.1 != factor.to.test],
                            factor.to.test)
      design.1 <- data.frame(design.1[,covariate.name.1])
    }else if(length(factor.to.test) == 3){
      # In this situation, user wants to compare two levels
      # make the second level as reference for convenience

      covariate.name.1 <- c(covariate.name.1[covariate.name.1 != factor.to.test[1]],
                            factor.to.test[1])
      design.1 <- data.frame(design.1[,covariate.name.1])
      # set the third element in 'factor.to.test' as reference level
      design.1[,factor.to.test[1]] <- relevel(design.1[,factor.to.test[1]], factor.to.test[3])
    }

  }

  # in our model, we also include cell type composition as main terms
  # in addition, design.2 provide covariates have NO cell type specific effects
  # so they are also modeled as main terms.
  # Put them in front of cell type specific covariates for convenience
  if(is.null(design.2)){
    dd <- design.1
  }else{
    dd <- cbind(design.2, design.1)
  }

  # create formula for main terms
  formul <- paste("~", paste(c(covariate.name.2, covariate.name.1), collapse = "+"))

  design_matrix <- model.matrix(as.formula(formul), dd)

  # remove replicated column
  design_matrix <- unique(design_matrix, MARGIN = 2)

  formul <- paste0("~ ", paste(colnames(design_matrix), sep = "+", collapse = "+"))
  return(list(design_matrix = design_matrix, design.1 = design.1, design.2 = design.2,
              all_coefs = c(covariate.name.1, covariate.name.2), formula = formul))
}
