# High-Dimensional Inference based on the Leave-One-Covariate-Out LASSO path
## Cao, X, Gregory, K.B, Wang, D

This R package is for high-dimensional inference using LOCO regularization path. Allows calculating **variable/feature importance** for regularized linear/logistic/GLM model , **varaible screening** for high-dimensional regression and **statistical inference** in the high-dimensional settings. 

**Parallel computing** is surpported to accelerate the inference process. 

The statistical inference procedure now are for high-dimensional linear regression only. Soon, we will add tools for GLM, Cox PH model and Gaussian graphical models.  

### Install the R package

To install, please use 
```R
devtools::install_github("devcao/LOCOpath")
```

### How to use?
#### Simple example
```R
require(LOCOpath)  ## load the package
set.seed(10)
data = depenDesign(n = 100, p = 12, beta = c(0,rep(1,2),rep(0,9)), rho = 0) ## generate some data
## calculate our test statistic for beta_1, using L1 norm
TS = ExactPath.TS(X=data$X, Y=data$Y, which.covariate = 1, betaNull = 0, multiTest = FALSE, norm = 'L1')
## calculate bootstrapped p-value for beta_1, using L1 norm 
rsmpl_rslt = Path.Resample(X=data$X, Y=data$Y, which.covariate = 1, betaNull = 0, multiTest = FALSE, norm = 'L1', B = 500, parallel = TRUE, beta.init = 'adaptive', beta.true = 0)
print(rsmpl_rslt$pval)
```
You should see: 
```R
> print(rsmpl_rslt$pval)
[1] 0.226
```

#### Power curve simulation
