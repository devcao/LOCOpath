#' Simulate the power of LOCO regularization path test statistics
#' @description This function simulate the power of LOCO regularization path test statistics
#' @param n: integar, sample size
#' @param p: integar, number of covariates
#' @param rho: correlation, if real value, $eqn{\Sigma_{i,j} = \rho ^ {|i-j|}}
#' 				if string, c("equl", "weak_equl")
#' @param beta: true beta coefficient
#' @param iter: number of iterations
#' @param B: number of bootstraps
#' @param which.covariate: args passed to Path.TS
#' @param betaNull: args passed to Path.TS
#' @param multiTest: args passed to Path.TS
#' @param ...: args passed to Path.TS
#' @return a vector of power for each Hypothesis test specified
#' @export
Path.Resample.Power = function(n = 100, p = 1000, beta=c(rep(1,10),rep(0,990)), rho=0.5, iter = 500, B = 500, setting = 'dep', which.covariate = 1, betaNull = 1, multiTest = FALSE, ...){
#Return the power fo Path-based test under different settings    
# Args:
#	setting: different settings, check 'pathwise_simu_setting.R for details'	
# 	rho: related to dependent design setting
#	n,p,beta : sample size, features, coefficients
# 	iter : # of iterations 
#
# Return:
#	Simulated power
#

	#TS = matrix(NA, iter, len(which.covariate))
	#b_size_la <-  matrix(0,iter,4)
  
  	path.power = array(0,dim = c(iter,length(which.covariate),4))

	for(s in 1:iter){
    
  		data = dataGen(setting = setting, n = n, p = p, beta = beta, rho = rho)

		X_sp = data$X
		Y_sp = data$Y
		
		results = Path.Resample(X = X_sp, Y = Y_sp, which.covariate = which.covariate, betaNull = betaNull, multiTest = multiTest, B = B, beta.true = beta, ...)

		######## keep track of MEMORY
		print("After Bootstrap:")
		print(mem_used())  
		########
		
		path.power[s,,] = results$rej

		if(s %% 10 == 0){  cat("Now Computing:",s,"\n") }
        
   
	}    

        
	path.power = apply(path.power,c(2,3),mean)
	
	return(path.power = path.power)
    
   
}    
    

