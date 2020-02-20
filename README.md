# A generalized framework for high-dimensional inference using Leave-One-Covariate-Out regularization path
## Cao, X, Gregory, K.B, Wang, D

#### For the **python** package, please go to [this link](https://github.com/statcao/LOCO_parse).

#### For the source code (includes both Python/C++/R code), please go to [this link](https://github.com/statcao/High-dimensional-inference).


This R package is for high-dimensional inference using LOCO regularization path. Allows calculating **variable/feature importance** for regularized linear/logistic/GLM model , **varaible screening** for high-dimensional regression and **statistical inference** in the high-dimensional settings. 

**Parallel computing** is surpported to accelerate the inference process. 


Our method is adaptable for all regularization techique, as long as you can solve the whole regularization path. We also provide the framework to feed in your own solution path. But keep in mind you need to verify the statistical correctness. (Our first paper verified it for linear models only). 

The statistical inference procedure now are for high-dimensional linear regression only. Soon, we will add tools for GLM, Cox PH model and Gaussian graphical models.  

If you need p-values for all variables, the current computing speed is slow, especily in the high-dimensional case. The core calculation is implemented in C++. Because we need to fit the whole path many times. the bottleneck of speed is the coordinate descent and Lars algorithm, which we used standard package glmnet, Lars and scikit learn.

### Install the R package

To install, please use 
```
devtools::install_github("devcao/LOCOpath")
```

### Simulation 
Check ```power_simulation_template.R``` for the code.



