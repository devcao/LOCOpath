


#' Generate data with design specified
#' @description This function wraps up different design for the simulated data 
#' @param setting: c("dep")
#' @param ...: args passed to depenDesign
#' @return a list with X and Y 
#' @export
dataGen = function(setting = 'dep', ...){
	if(setting == 'dep'){

		return(depenDesign(...))

		}else if (setting == 'other'){
		continue 
		#### add more settings here
	}

}


#' Generate data with correlated design 
#' @description This function Generate the data for simulation
#' @param n: integar, sample size
#' @param p: integar, number of covariates
#' @param rho: correlation, if real value, $eqn{\Sigma_{i,j} = \rho ^ {|i-j|}}
#' 				if string, c("equl", "weak_equl")
#' @return a list with X and Y
#' @export
depenDesign = function(n, p, beta, rho){

	if(rho == 'equl'){  # equi corr
    	Sigma = matrix(rep(0.8,p*p),p,p)
    	diag(Sigma) = rep(1,p)

    	Mu=rep(0,p)
		X=rmvnorm(n,mean=Mu,sigma=Sigma)
		Y <- X %*% beta + rnorm(n,0,1)   
   
		}else if (rho == "weak_equl"){
				
				Sigma = matrix(rep(0.5,p*p),p,p)
    			diag(Sigma) = rep(1,p)

    			Mu=rep(0,p)
				X=rmvnorm(n,mean=Mu,sigma=Sigma)
				Y <- X %*% beta + rnorm(n,0,1)   
   
			}else if (rho > 0){  # Topliez corr

    		Sigma = diag(p)

    		for(i in 1:p){
      			for(j in 1:p){
        			Sigma[i,j]=rho^(abs(i-j))
        		}
    		}

    		Mu = rep(0,p)
			X = rmvnorm(n,mean=Mu,sigma=Sigma)
			Y <- X %*% beta + rnorm(n,0,1)  

			}else if (rho == 0){  # independent

				X <- matrix(rnorm(n*p),nrow=n)
				Y <- X %*% beta + rnorm(n,0,1)
			
			}


	return(list(X = X, Y = Y))   

}


# Setting 2