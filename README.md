# High-Dimensional Inference based on the Leave-One-Covariate-Out LASSO path
## Cao, X, Gregory, K.B, Wang, D

This R package is for high-dimensional inference using LOCO regularization path. Allows calculating **variable/feature importance** for regularized linear/logistic/GLM model , **varaible screening** for high-dimensional regression and **statistical inference** in the high-dimensional settings. 

**Parallel computing** is surpported to accelerate the inference process. 

The statistical inference procedure now are for high-dimensional linear regression only. Soon, we will add tools for GLM, Cox PH model and Gaussian graphical models.  

### Install the R package

To install, please use 
```
devtools::install_github("devcao/LOCOpath")
```

### How to use?
#### Simple example
```
require(LOCOpath)
data = depenDesign(n = 100, p = 1000, beta = c(0,rep(1,9),rep(0,990)), rho = 0)

```
