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
set.seed(10) ## set seed
data = depenDesign(n = 100, p = 12, beta = c(0,rep(1,2),rep(0,9)), rho = 0) ## generate some data

## calculate our test statistic for beta_1, using L1 norm
TS = ExactPath.TS(X = data$X, Y = data$Y, which.covariate = 1, betaNull = 0, multiTest = FALSE, norm = 'L1')
## calculate bootstrapped p-value for beta_1, using L1 norm 
rsmpl_rslt = Path.Resample(X = data$X, Y = data$Y, 
                           which.covariate = 1, betaNull = 0, multiTest = FALSE, norm = 'L1', 
                           B = 500, parallel = TRUE, beta.init = 'adaptive', beta.true = 0)
print(rsmpl_rslt$pval)
```
You should the following p-value for beta_1
```R
> print(rsmpl_rslt$pval)
[1] 0.226
```

#### Power curve simulation 
```
require(LOCOpath)
## model setting 
n = 100; p = 12; rho = 0
## Power: 500 iterations; Bootstrap: 500 replicates
iter = 500; B = 500
power_indep_L1 = matrix(0, 11, 4)
i = 1
for (b in seq(0,1,0.1)){
  power_indep_L1[i, ] = Path.Resample.Power(n = n, p = p, beta=c(0+b,rep(1,2),rep(0,9)), 
                                     rho=rho, iter = iter, B = B, setting = 'dep', 
                                     which.covariate = 1, betaNull = 0, multiTest = FALSE, 
                                     parallel = TRUE, norm = 'L1', path.method ='lars', beta.init = 'adaptive')
}
### and also calculate the power curve with rho = c(0, 0.5, 0.9) and norm = c('L1', 'L2.squared', 'L_inf')
## plot power curve 
require(latex2exp)
cex = 1.4
par(mfrow=c(3,1),mar=c(0,0,0,0),oma = c(5.1, 4.1, 4.1, 2.1))


plot(NA,xlim=c(0,1),ylim=c(0,1),xaxt="n")
lines(power_indep_L1[,3]~beta, pch = 19, type = 'b', lwd = 2)
lines(power_indep_L_inf[,3]~beta, pch = 1, type = 'b', lwd = 2)
lines(power_indep_t_test[,3]~beta, pch = 2, type = 'b', lwd = 2)
abline(h = 0.05, lty = 2)

legend( x = grconvertX(0, from = 'nfc', to = 'user'),
        y = grconvertX(1, from = 'nfc', to = 'user'),
        legend=TeX('$\\Sigma = \\mathbf{I}_p$ '),
        bty="n", cex = cex)

plot(NA,xlim=c(0,1),ylim=c(0,1),xaxt="n")
lines(power_rho05_L1[,3]~beta, pch = 19, type = 'b', lwd = 2)
lines(power_rho05_L_inf[,3]~beta, pch = 1, type = 'b', lwd = 2)
lines(power_rho05_t_test[,3]~beta, pch = 2, type = 'b', lwd = 2)
abline(h = 0.05, lty = 2)

legend( x = grconvertX(0, from = 'nfc', to = 'user'),
        y = grconvertX(1, from = 'nfc', to = 'user'),
        legend=TeX('$\\Sigma = (0.5^{|i-j|})_{1\\leq i , j \\leq p}$ '),
        bty="n", cex = cex)

plot(NA,xlim=c(0,1),ylim=c(0,1),xaxt="n")
lines(power_rho09_L1[,3]~beta, pch = 19, type = 'b', lwd = 2)
lines(power_rho09_L_inf[,3]~beta, pch = 1, type = 'b', lwd = 2)
lines(power_rho09_t_test[,3]~beta, pch = 2, type = 'b', lwd = 2)
abline(h = 0.05, lty = 2)

legend( x = grconvertX(0, from = 'nfc', to = 'user'),
        y = grconvertX(1, from = 'nfc', to = 'user'),
        legend=TeX('$\\Sigma = (0.9^{|i-j|})_{1\\leq i , j \\leq p}$ '),
        bty="n", cex = cex)

mtext(side=1, outer=TRUE, expression(beta[1]),line=3)
mtext(side=2, outer=TRUE, "Empirical power",line=2.5)


legend( x = grconvertX(0.2, from = 'nfc', to = 'user'),
        y = grconvertX(3.1, from = 'nfc', to = 'user'),
        legend=c( expression(T[1](1,1)),
                  expression(T[1](infinity,infinity)),
                  "T-test"),
        pch=c(19,1,2),
        lwd = 1,
        bty="n", xpd = NA, horiz = TRUE, cex = cex-0.2)



```
![Optional Text](./new_12.pdf)
