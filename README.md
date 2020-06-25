# LOCOpath
This R package is for high-dimensional inference using Leave-One-Covariate-Out(LOCO) LASSO path.

This R package accompanies the paper:

Cao, X., Gregory, K. B., & Wang, D. High-Dimensional Inference based on the Leave-One-Covariate-Out LASSO path

## Install the R package

To install, please first install R package ```devtools``` and then 
```R
devtools::install_github("devcao/LOCOpath")
```
#### Some recent issue in installing
**R removed one of our dependence R package ```parcor```**. To install, you may need to manually install ```parcor``` first.
To install ```parcor```, first install ```parcor``` dependence package
```R
install.packages("Epi")
install.packages("GeneNet")
```
Then download ```ppls``` and ```parcor```. And install 
```R
install.packages('YOUR_PATH/ppls_1.6-1.1.tar.gz', repos = NULL, type="source")
install.packages('YOUR_PATH/parcor_0.2-6.tar.gz', repos = NULL, type="source")
```
Now you should be able to install our packages!

## How to use
Feel free to check out [some working examples](./Examples.md) for more detailed examples.
