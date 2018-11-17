


#' The engine to perform residual bootsrrap of the Leave-One-Covariate-Out regularization path test staitsic
#' @description This function performs residual bootsrrap the Leave-One-Covariate-Out regularization path test staitsic
#' @param X: a n by p matrix, design matrix 
#' @param Y: a n by 1 matrx, response 
#' @param multiTest: boolean, TRUE: test simultaneously, FALSE: test individually 
#' @param betaNull: same args for ExactPath.TS
#' @param B: integer, number of bootstraps 
#' @param parallel: boolean, run the whole bootstrap process in parallell, useful if run in a multicore cpu, especially in cluster
#' @param exact: boolean, using ExactPath.TS or AprroxPath.TS
#' @param b.Null: Args from Path.Resample
#' @param beta.index: Args from Path.Resample
#' @param ...: args wull be passed to ExactPath.TS or ApproxPath.Ts
#' @return Will be passed to Path.Resample
#' @export


Path.Resample.Process = function(X, Y, multiTest, residual, b.Null, beta.index, betaNull, B = 500, exact, parallel = FALSE, ...){
# Bootstrap the null distribution of Path-based statistic of coef beta.index
#
# Args:
#	X, Y, which.covariate : feed in to Path-based TS function
#	B : # of bootstrap replications
#	parallel : run in parallel or nor
#	exact : use exact TS or approx TS
#	beta.index : which coef 
#
# Return:  
#	A vector of the bootstrapped null 
	n = nrow(X)
	p = ncol(X)

	TS_null = numeric()
  
  	if(parallel){ # running in parallel
    	mat = list()


    	for(bs in 1:B){
      		ind = sample(1:n,replace = TRUE)
      		boot_residual = residual[ind]
      		#b_null = bhat
      		#b_null[beta.index] = 0
      		Y = X %*% b.Null + boot_residual
      		mat[[bs]] = cbind(X,Y) 
    	}
  		
  		#rgs = list(...)
    	#Args = c(which.covariate = beta.index, betaNull = betaNull, exact = exact, multiTest = multiTest, args)

    	# On a cluster, just use
    	
    	no_cores <- detectCores() 
    	cat("n_cores detected:", no_cores, "\n")
    	# Initiate cluster
    	#cl <- makeCluster(no_cores)
    	cl <- makeCluster(no_cores, type = "FORK")
    	###### if using window ###### 
      # load special packages
    	#clusterEvalQ(cl, .libPaths("~/R"))
    	#clusterEvalQ(cl, library(glmnet))
    	#clusterEvalQ(cl, library(lars))
    	#clusterEvalQ(cl, library(MASS))
    	#clusterEvalQ(cl,library(pryr))
    	#clusterEvalQ(cl,library(plus))
    	#clusterEvalQ(cl,source("~/hdi_path/bin/pathwise_ts.R"))

    	#clusterExport(cl, varlist = c("beta.index", "exact", "betaNull", "multiTest",...), envir = environment())
      	#clusterExport(cl, varlist = 'Args', envir = environment())
      ###### if using window ######

    	re_list = parLapply(cl, mat, Path.TS.Para, exact = exact, multiTest = multiTest, which.covariate = beta.index, betaNull = betaNull, ...)
    	#re_list = parLapply(cl, mat, Path.TS.Para, list = Args)
    
    	######## in case run out of MEMORY
    	print("Cluster MEM:")
    	print(mem_used())
    	########
    	stopCluster(cl)  # END parallel bootsrap


    	for(bss in 1:B){

      		TS_null[bss] = re_list[[bss]]
      
    	}

    	return(TS_null)

    }else{ # not parallel, could be slow

  		for(bs in 1:B){   
    		ind = sample(1:n,replace = TRUE)
    		boot_residual = residual[ind]
    		#b_null = bhat 
    		#b_null[beta.index] = 0
    		Y = X %*% b.Null + boot_residual
    		
    		TS_null[bs] = Path.TS(exact = exact, X = X, Y = Y, multiTest = multiTest, which.covariate = beta.index, betaNull = betaNull,...)
	     	
  		}

  		return(TS_null)
  	}
	
    

}

    






	

    
    
   

    	
	
