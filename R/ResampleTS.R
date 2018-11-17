

#' Residual bootsrrap the Leave-One-Covariate-Out regularization path test staitsic
#' @description This function performs residual bootsrrap the Leave-One-Covariate-Out regularization path test staitsic and would return test p-value
#' @param X: a n by p matrix, design matrix 
#' @param Y: a n by 1 matrx, response 
#' @param which.covariate: a vector or a list of vector, specify which covaritate \eqn{ \beta_j } to perform our test 
#' @param betaNull: a vector or a list of vector, should be the same shape as which.covariate. 
#'					Specify the Null hypothesis, \eqn{H_0: \beta_j = betaNull vs H_1 \neq betaNull }
#' 					if using multiTest = FALSE, then will return 				
#' @param multiTest: boolean, TRUE: test simultaneously, FALSE: test individually 
#' @param B: integer, number of bootstraps 
#' @param parallel: boolean, run the whole bootstrap process in parallell, useful if run in a multicore cpu, especially in cluster
#' @param exact: boolean, using ExactPath.TS or AprroxPath.TS
#' @param beta.init: str, c("adaptive", "de-sparse", "MC+", "SCAD", "Truth"). Under developing, please use "adaptive" and "Truth". 
#'						Needed if perform high-dimensional regression, low-dimensional regression will use LSE.
#' @param beta.true: true value of beta, if your beta.init = "Truth", need to specify this
#' @param ...: args wull be passed to ExactPath.TS or ApproxPath.Ts
#' @return a list of p-value, TS, TS_null, rej or nor
#' @examples
#' X = matrix(rnorm(100*12), 100, 12)
#' beta = c(1,2, rep(0, 10))
#' Y = X %*% beta + rnorm(100)
#' # Test \eqn{H_0: \beta_1 = 1} and \eqn{H_0: \beta_1 = 2} individually 
#' Path.Resample(X, Y, which.covariate = c(1, 2), betaNull = c(1, 2), multiTest = FALSE, B = 100)
#' # Test \eqn{H_0: \beta_1 = 1, \beta_1 = 2} simultaneously
#' Path.Resample(X, Y, which.covariate = list(c(1, 2)), betaNull = list(c(1, 2)), multiTest = TRUE, B = 100)
#' # Test \eqn{H_0: \beta_1 = 1, \beta_2 = 2} simultaneously and test \eqn{H_0: \beta_1 = 0, \beta_2 = 0}
#' Path.Resample(X, Y, which.covariate = list(c(1, 2), c(1, 2)), betaNull = list(c(1, 2), c(0, 0)), multiTest = TRUE, B = 100)
#' @export



Path.Resample = function(X, Y, which.covariate, betaNull, multiTest, B = 500, parallel = FALSE, exact = TRUE, beta.init = 'adaptive', beta.true = beta, ...){
# Bootstrap the null distribution of Path-based statistic, and return reject or not
#
# Args:
#	X, Y, which.covariate : feed in to Path-based TS function
#	B : # of bootstrap replications
#	parallel : run in parallel or nor
#	exact : use exact TS or approx TS
#	beta.true : for simulation only, the true value of beta would be used as initial estimates.
#	
# Return:  
#	Reject or not under alpha = 0.2,0.1,0.05,0.01
# 	p.values of the test
  
	
  n = nrow(X)
  p = ncol(X)
  rej = matrix(0,length(which.covariate),4)  # 1st : which cov, 2nd: na, 3rd: which alpha, 0.2,0.1,0.05,0.01
  pval = numeric()

  TS = Path.TS(exact = exact, X = X, Y = Y, which.covariate = which.covariate, betaNull = betaNull, multiTest = multiTest,...)
  
  



  if(p >= n){ # high dimension we can try...

  	if(beta.init == "adaptive"){
  		bhat = adalasso(X = X, y = Y, k = 10, use.Gram = FALSE,both = TRUE, intercept = FALSE)$coefficients.adalasso
  
    	}else if (beta.init == "de-sparse"){
    		
    		bhat = as.vector(lasso.proj(X, Y, standardize = TRUE, parallel = TRUE, ncores = 40)$bhat)

    		}else if (beta.init == "MC+"){
    			
    			bhat = coef(cv.ncvreg(X = X, y = Y, penalty = "MCP",family = "gaussian", nfold= 10))[-1]
          
    			}else if (beta.init == "SCAD"){
    				bhat = coef(cv.ncvreg(X = X, y = Y, penalty = "SCAD",family = "gaussian", nfold= 10))[-1]

    			}else if (beta.init == "Truth"){
    				bhat = beta.true
    			} 

    
    residual = Y - X%*%bhat
    
  	}else{ # low dimenstion just use LSE
    	bhat = ginv(t(X)%*%X)%*%t(X)%*%Y    
    	residual = Y - X%*%bhat
  	
  	}
  
  #TS_null = matrix(NA, nrow = B, ncol = length(which.covariate))

  ################### HERE WE GO ! ! ! ###########################################
  	

  	###################### This part could be parallelized ##################	
  	count = 1

  	for(wc_cov in which.covariate){   

  		b.Null = bhat
  		#b.Null[wc_cov] = 0
 
  		if(multiTest) { 
  			to.which.covariate = list(wc_cov)
  			to.betaNull = list(betaNull[[count]])

  			b.Null[wc_cov] = betaNull[[count]]


  		}else{
  			to.which.covariate = wc_cov
  			to.betaNull = betaNull[count]

  			b.Null[wc_cov] = betaNull[count]

  		} # then run multiple testing


  		TS_null = Path.Resample.Process(X = X, Y = Y, multiTest = multiTest, residual = residual, b.Null = b.Null, betaNull = to.betaNull, 
  												beta.index = to.which.covariate, B = B, exact = exact, parallel = parallel, ...)

  		rej[count,1] = TS[count] > quantile(TS_null,0.8)
    	rej[count,2] = TS[count] > quantile(TS_null,0.9)
    	rej[count,3] = TS[count] > quantile(TS_null,0.95)
    	rej[count,4] = TS[count] > quantile(TS_null,0.99)
    	pval[count] = mean(TS_null > TS[count])

    	count = count + 1
    	
  	}


  	
  
  ##########################################################

  return(list(rej = rej, pval = pval, TS_null = TS_null, TS = TS))
  
} 

