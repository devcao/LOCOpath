
#' Approximate Calculation of Leave-One-Covariate-Out regularization path test staitsic
#' @description This function performs the approximate Leave-One-Covariate-Out regularization path estimator using interpolation and approximate integration
#' @param X: a n by p matrix, design matrix 
#' @param Y: a n by 1 matrx, response 
#' @param which.covariate: a vector or a list of vector, specify which covaritate \eqn{ \beta_j } to perform our test 
#' @param betaNull: a vector or a list of vector, should be the same shape as which.covariate. 
#'					Specify the Null hypothesis, \eqn{H_0: \beta_j = betaNull vs H_1 \neq betaNull }
#' 					if using multiTest = FALSE, then will return 				
#' @param multiTest: boolean, TRUE: test simultaneously, FALSE: test individually 
#' @param path.method: string, c("lars", "plus.lasso", "plus.mc+", "plus.scad" )
#' @param norm: string, c("L1", "L2", "L2.squared", "L_inf")
#' @param normalize: boolean, using normalization or not
#' @param intercept: boolean, include intercept or not
#' @return a vector of test statistic for each test specified
#' @examples
#' X = matrix(rnorm(100*12), 100, 12)
#' beta = c(1,2, rep(0, 10))
#' Y = X %*% beta + rnorm(100)
#' # Test \eqn{H_0: \beta_1 = 1} and \eqn{H_0: \beta_1 = 2} individually 
#' ApproxPath.TS(X, Y, which.covariate = c(1, 2), betaNull = c(1, 2), multiTest = FALSE)
#' # Test \eqn{H_0: \beta_1 = 1, \beta_1 = 2} simultaneously
#' ApproxPath.TS(X, Y, which.covariate = list(c(1, 2)), betaNull = list(c(1, 2)), multiTest = TRUE)
#' # Test \eqn{H_0: \beta_1 = 1, \beta_2 = 2} simultaneously and test \eqn{H_0: \beta_1 = 0, \beta_2 = 0}
#' ApproxPath.TS(X, Y, which.covariate = list(c(1, 2), c(1, 2)), betaNull = list(c(1, 2), c(0, 0)), multiTest = TRUE)
#' @export
ApproxPath.TS <- function(X,Y, which.covariate, betaNull, multiTest, path.method = "lars", norm = 'L2.squared',normalize = TRUE, intercept = FALSE, ridgePen = 0, numInterpolation = 10000){
# Calculate PATH statistic using interpolating approximation 
#
# Args:
#	X,Y: design matrix and response vector
#	which.covariate: a vector, indicating which covariate we will be computing
# 	normalize: argguments of lars 
# 	nrom: indludes "L1", "L2", "L2.squared","L_inf"
# 	ridgePen: ridge penalty term, 0 indicates no ridge penalty
# 	numInterpolation: number of interpolations, default 10000.
#betaNull: same size and same data type with which.covariate, specify the null hypothesis H0: beta = betaNull. 'Zero' or user sepcify,
#		like betaNull = c(1,1,1), which.covariate = c(1,2,3),betaNull = list(c(1,1,1),c(0,0)), which.covariate = list(c(1,2,3) ,c(5,6))

# Returns:
#	A vector of PATH statistic
#
#	

	n = nrow(X)
	p = ncol(X)
	
	if (ridgePen > 0){
		zero_vec = rep(0,p)
		Iden_mat = diag(rep(1,p))
		X = rbind(X,ridgePen*Iden_mat) 
		Y = append(Y,zero_vec)
	
	}
################# Here we go ! ! ! ##############################################



	l <- 1
	TS <- numeric(length(which.covariate))


	for (j in which.covariate){
 
		# first condition : whether test all 0 or nor
 		
		if(path.method == "lars"){
				# 2nd condition, muti-test or not 
			if(multiTest & is.list(which.covariate) & is.list(betaNull)){
					

					adjust.X = rowSums(  t( apply(X[,j], 1, function(x){return(x*betaNull[[l]])}) )  )
					newY = Y - adjust.X

					X.sc = scale(X)


					lars.out <- lars(X.sc, newY, type = "lasso", intercept = intercept, use.Gram = FALSE, normalize=normalize)
					lambda.hat <- sort(c(lars.out$lambda,0),decreasing=FALSE)
					beta.hat <- coef(lars.out)[seq(length(lambda.hat),1,-1),] 

					lars.j.out <- lars(X.sc[, -j], newY, type = "lasso", intercept = intercept, use.Gram = FALSE,normalize=normalize)
					lambda.j.hat <- sort(c(lars.j.out$lambda,0),decreasing=FALSE)
					beta_val <- coef(lars.j.out)[seq(length(lambda.j.hat),1,-1),]

			
				}else if((!multiTest)){  #indivdual test
						
					newY = Y - betaNull[l] * X[,j]	
					X.sc =scale(X)

					lars.out <- lars(X.sc, newY , type = "lasso", intercept = intercept, use.Gram = FALSE, normalize=normalize)
					lambda.hat <- sort(c(lars.out$lambda,0),decreasing=FALSE)
					beta.hat <- coef(lars.out)[seq(length(lambda.hat),1,-1),] 

					lars.j.out <- lars(X.sc[, -j], newY, type = "lasso", intercept = intercept, use.Gram = FALSE,normalize=normalize)
					lambda.j.hat <- sort(c(lars.j.out$lambda,0),decreasing=FALSE)
					beta_val <- coef(lars.j.out)[seq(length(lambda.j.hat),1,-1),]

				}else{
					stop("wrong input")
				}	

			}else if(path.method == "plus.lasso"){


				if(multiTest & is.list(which.covariate) & is.list(betaNull)){
					

					adjust.X = rowSums(  t( apply(X[,j], 1, function(x){return(x*betaNull[[l]])}) )  )
					newY = Y - adjust.X

					X.sc = scale(X)

					#### testing
					lars.out <- plus(X.sc, as.vector(newY), method = "lasso",intercept = intercept, use.Gram = FALSE, normalize=normalize)
					lambda.hat <- sort(lars.out$lam.path, decreasing=FALSE)
					beta.hat <- coef(lars.out, lam = lambda.hat)

					lars.j.out <- plus(X.sc[, -j], as.vector(newY), method = "lasso", intercept = intercept, use.Gram = FALSE,normalize=normalize)
					lambda.j.hat <- sort(lars.j.out$lam.path, decreasing=FALSE)
					beta_val <- coef(lars.j.out, lam = lambda.j.hat)

			
				}else if((!multiTest)){  #indivdual test
						
					newY = Y - betaNull[l] * X[,j]	
					X.sc =scale(X)


					lars.out <- plus(X.sc, as.vector(newY) , method = "lasso", intercept = intercept, use.Gram = FALSE, normalize=normalize)
					lambda.hat <- sort(lars.out$lam.path, decreasing=FALSE)
					beta.hat <- coef(lars.out, lam = lambda.hat)

					

					lars.j.out <- plus(X.sc[, -j], as.vector(newY), method = "lasso", intercept = intercept, use.Gram = FALSE,normalize=normalize)
					lambda.j.hat <- sort(lars.j.out$lam.path, decreasing=FALSE)
					beta_val <- coef(lars.j.out, lam = lambda.j.hat)

				}else{
					stop("wrong input")
				}	

				
			}else if (path.method == "plus.mc+"){
				if(multiTest & is.list(which.covariate) & is.list(betaNull)){
					

					adjust.X = rowSums(  t( apply(X[,j], 1, function(x){return(x*betaNull[[l]])}) )  )
					newY = Y - adjust.X

					X.sc = scale(X)


					lars.out <- plus(X.sc, as.vector(newY), method = "mc+", intercept = intercept, use.Gram = FALSE, normalize=normalize)
					lambda.hat <- sort(lars.out$lam.path,decreasing=FALSE)
					beta.hat <- coef(lars.out, lam = lambda.hat)
					
					lars.j.out <- plus(X.sc[, -j], as.vector(newY), method = "mc+", intercept = intercept, use.Gram = FALSE,normalize=normalize)
					lambda.j.hat <- sort(lars.j.out$lam.path,decreasing=FALSE)
					beta_val <- coef(lars.j.out, lam = lambda.j.hat)
				
				}else if((!multiTest)){  #indivdual test
						
					newY = Y - betaNull[l] * X[,j]	
					X.sc =scale(X)

					lars.out <- plus(X.sc, as.vector(newY) , method = "mc+", intercept = intercept, use.Gram = FALSE, normalize=normalize)
					lambda.hat <- sort(lars.out$lam.path,decreasing=FALSE)
					beta.hat <- coef(lars.out, lam = lambda.hat)

					lars.j.out <- plus(X.sc[, -j], as.vector(newY), method = "mc+", intercept = intercept, use.Gram = FALSE,normalize=normalize)
					lambda.j.hat <- sort(lars.j.out$lam.path,decreasing=FALSE)
					beta_val <- coef(lars.j.out, lam = lambda.j.hat)

				}else{
					stop("wrong input")
				}




			}else if (path.method == "plus.scad"){

				if(multiTest & is.list(which.covariate) & is.list(betaNull)){
					

					adjust.X = rowSums(  t( apply(X[,j], 1, function(x){return(x*betaNull[[l]])}) )  )
					newY = Y - adjust.X

					X.sc = scale(X)


					lars.out <- plus(X.sc, as.vector(newY), method = "scad", intercept = intercept, use.Gram = FALSE, normalize=normalize)
					lambda.hat <- sort(lars.out$lam.path,decreasing=FALSE)
					beta.hat <- coef(lars.out, lam = lambda.hat)

					lars.j.out <- plus(X.sc[, -j], as.vector(newY), method = "scad", intercept = intercept, use.Gram = FALSE,normalize=normalize)
					lambda.j.hat <- sort(lars.j.out$lam.path,decreasing=FALSE)
					beta_val <- coef(lars.j.out, lam = lambda.j.hat)

					
					
			
				}else if((!multiTest)){  #indivdual test
						
					newY = Y - betaNull[l] * X[,j]	
					X.sc =scale(X)

					lars.out <- plus(X.sc, as.vector(newY) , method = "scad", intercept = intercept, use.Gram = FALSE, normalize=normalize)
					lambda.hat <- sort(lars.out$lam.path,decreasing=FALSE)
					beta.hat <- coef(lars.out, lam = lambda.hat)

					lars.j.out <- plus(X.sc[, -j], as.vector(newY), method = "scad", intercept = intercept, use.Gram = FALSE,normalize=normalize)
					lambda.j.hat <- sort(lars.j.out$lam.path,decreasing=FALSE)
					beta_val <- coef(lars.j.out, lam = lambda.j.hat)
				
				}else{
					stop("wrong input")
				}




			}


			#### remove the unoverlapping part (test)

		if(lambda.hat[1] != lambda.j.hat[1]){

			leftmost = c(lambda.hat[1], lambda.j.hat[1])
			whichone = which.max(leftmost)

			if(whichone ==1){

				lambda.j.hat <- lambda.j.hat[lambda.j.hat >= lambda.hat[1]]
				lambda.j.hat <- c(lambda.hat[1], lambda.j.hat)
			}else{

				lambda.hat <- lambda.hat[lambda.hat >= lambda.j.hat[1]]
				lambda.hat <- c(lambda.j.hat[1], lambda.hat)
			
			}
			
			beta.hat <- coef(lars.out, lam = lambda.hat)
			beta_val <- coef(lars.j.out, lam = lambda.j.hat)
					
			
		}

		####


		new_beta <- matrix(0, dim(beta_val)[1], p)
		new_beta[, -j] <- beta_val
		beta.j.hat <- new_beta
 
		union.lambda <- sort(unique(c(lambda.hat,lambda.j.hat)),decreasing=FALSE)  # union of lambda values from both solution paths.
		M <- length(union.lambda) # cardinality of union of lambda values from both solution paths
 
		beta.hat.union.lambda <- beta.j.hat.union.lambda <- matrix(NA,length(union.lambda),p)
 
		TS.k <- numeric()
		
 		for (k in 1:p){

 			# using logrithmic scale: 
			#interpolating.lambda = exp(seq(log(min(union.lambda[union.lambda > 0])),log(max(union.lambda)),length.out = numInterpolation+1))   
			# using linear scale: 
			interpolating.lambda = seq(min(union.lambda),max(union.lambda),length.out = numInterpolation+1)   
		
			approx.beta.hat = approx(x = lambda.hat, y = beta.hat[,k], xout = interpolating.lambda,yright = 0)
			approx.beta.j.hat = approx(x = lambda.j.hat, y = beta.j.hat[,k], xout = interpolating.lambda,yright = 0)

			#beta.hat.union.lambda[,k] <- approx.beta.hat$y 
			#beta.j.hat.union.lambda[,k] <- approx.beta.j.hat$y
			# get absolute difference between beta.hat and beta.j.hat at all values of lambda in union lambda
			
			if (norm == "L2.squared"){
				#delta <- (beta.hat.union.lambda[,k] - beta.j.hat.union.lambda[,k])
				#TS.k[k] <- sum(diff(union.lambda)*(delta[-M]^2 + diff(delta) * delta[-M] + (1/3) * diff(delta)^2 ))
				
				TS.k[k] <- sum(((approx.beta.hat$y - approx.beta.j.hat$y)^2)[-1] * diff(interpolating.lambda))

			}else if (norm == "L1" ){
				TS.k[k] <- sum((abs(approx.beta.hat$y - approx.beta.j.hat$y))[-1] * diff(interpolating.lambda))

			}else if (norm == "L2"){
				TS.k[k] <- sum(((approx.beta.hat$y - approx.beta.j.hat$y)^2)[-1] * diff(interpolating.lambda))

			}else if (norm == "L_inf"){
				TS.k[k] <- max(abs(approx.beta.hat$y - approx.beta.j.hat$y))
				

			}
			

		}

		

		if (norm == "L2.squared" | norm == "L1"){
			TS[l] <- sum(TS.k)	
		}else if (norm == "L2"){
			TS[l] = sqrt(sum(TS.k))
		}else if (norm == "L_inf"){
			TS[l] <-  max(TS.k)
		}

 
		
 
		l <- l + 1





 
	}
 
	return(TS)
 
}

