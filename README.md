# A unified framework for high-dimensional inference using LOCO regularization path
## Cao, X, Gregory, K.B, Wang, D
### Maintained by Xiangyang Cao

This R package is for high-dimensional inference using LOCO regularization path, contains tools for computing test statistics and p-values. Parallel computing is surpported to accelerate the inference process. 

Our method is adaptable for all regularization techique, as long as you can solve the whole regularization path. We also provide the framework to feed in your own solution path. But keep in mind you need to verify the statistical correctness. (Our first paper verified it for linear models only). 

We are now focus on the inference for high-dimensional linear regression. Soon, we will add tools for GLM, Cox PH model and Gaussian graphical models.  

If you need p-values for all variables, the current computing speed is slow, especily for high-dimensional case. We are optimizing the performance of this package using Rcpp. To conduct multiple tests, we recommend use our simultaneous test.

### Install the R package

To install, please use 
```
install_github("statcao/LOCOpath")
```

