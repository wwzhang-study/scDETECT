# scDETECT: single-cell RNA-seq DE analysis accounting for the cell type correlations

`scDETECT` is an R package for single-cell RNA-seq DE analysis accounting for the cell type correlations. 
It first constructs a cell type hierarchical tree to represent the cell type correlations, and then incorporates the correlation in the calculation of the prior probabilities for DE.
After parameter estimation in data model, the posterior probabilities are then derived 
based on the Bayesian model for determining the DE state of each gene.
The main intuition of the method is that the inference for DE states for a gene in one 
cell type would be influenced by other cell types. Jointly modeling the expression 
in all cell types creates information sharing, which will improve the statistical inference.


## Installation
You can install `scDETECT` from [GitHub](https://github.com/wwzhang-study/scDETECT) using the [`devtools`](https://cran.r-project.org/web/packages/devtools/index.html) package:

```R
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")
    
devtools::install_github('wwzhang-study/scDETECT', dependencies=T, build_vignettes = T)

library(scDETECT)
```

## Getting started
In the `scDETECT` package, we provide a toy simulated dataset that you can directly 
generate using the following code:

```R
N <- 40
K <- 3
C <- 10
P <- 500
design <- data.frame(disease=factor(sample(0:1,size = N,replace=TRUE)))
Y = list()
res = list()
for (i in 1:K){
  y.pseudo = c()
  for (j in 1:N){
    y = matrix(rnbinom(P*C, size = 1, mu = 1), nrow = P, byrow = FALSE)
    y.pseudo = cbind(y.pseudo, rowSums(y))
  }
  Y[[i]] = y.pseudo
  rownames(Y[[i]]) <- paste0('gene',1:P)
}
```

In the simulation, we generate 3 cell types and 20 samples per group. In each cell type-sample combination,
we simulate single cell RNA-seq data for 500 genes and 10 cells for each cell type of each sample. Then we get summation to get pseudo-bulk data `Y` for each cell type.  

Next, you can run `scDETECT` with a single line of code:

```R
res_scDETECT <- scDETECT(Y_raw = Y, design.1 = design,design.2 = NULL,factor.to.test = 'disease',cutoff.tree = c('tstat',2.58),cutoff.prior.prob = c('pval',0.01))
```

Some explanations about the parameters:

- **design.1` and `design.2:** Covariates representing interested factors to be tested.
- **factor.to.test:** A phenotype name, e.g. "disease". 

- **cutoff.tree:** A cut off used to define DE state to estimate tree, which could be `fdr`, `pval` or `stat`.

- **cutoff.prior.prob:** A cut off used to define DE state to estimate prior probability of nodes on tree, which could be `fdr` or `pval`. 

- **pval:** A matrix of p-values from DESeq2. 

- **p.adj:** A matrix of adjusted p-values from DESeq2.

- **tree:** A hierarchical tree structure used to account cell type correlation.

- **p.matrix.input:** A matrix of prior probability on each node of the tree structure.

- **de.state:** DE state of each feature in each cell type.

- **similarity.function:** A custom function used to calculate similarity between cell types that used for tree structure estimation.

- **parallel.core:** The number of cores for parallel running.

- **corr.fig=TRUE:** Correlation between cell types will be plotted using function `plotCorr()`. 

- **run.time=TRUE:** Running time in seconds will be reported.

- **tree.type:** The tree type for inference.

We can have posterior probability of DE for each gene in each cell type:

```R
head(res_scDETECT$tree_res$full$pp)
```

The correlation between cell types was captured by a hierarchical tree estimated
from test statistics of DESeq2 result:
```R
res_scDETECT$tree_res$full$tree_structure
```

Cell types with smaller distance means they are stronger correlated. Different tree structures could be customized. Another simpler tree structure is also used for inference:

```R
res_scDETECT$tree_res$single$tree_structure

For more details about how to run `scDETECT` and classic DE methods, please refer to `vignette("scDETECT")`.


## Citation
Coming soon.
