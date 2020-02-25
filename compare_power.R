##################################################################################

# # Includes Power simulation of Projection method, etc.

##################################################################################




#######
##### loading required packages
#setwd('/hdi_simu')
require(LOCOpath)
len = length
######



################ de-sparsified lasso #########################
desparse.Power = function(n = 100, p = 1000, beta, rho, iter = 500, setting = 'dep', which.covariate, betaNull){
#Return the power of de-sparsified under different settings    
# Args:
#	setting: different settings, check 'pathwise_simu_setting.R for details'	
# 	rho: related to dependent design setting
#	n,p,beta : sample size, features, coefficients
# 	iter : # of iterations 
#
# Return:
#	A matrix of len(which.covariate) * 4 : Simulated power for j-th covariate
#  
  
  if(len(betaNull) > 1){stop("now only support compute power of 1 coefficients")}

  proj.power = matrix(0,len(which.covariate),4)
  
  pval = matrix(NA, iter, p)
  
  for(s in 1:iter){
    data = dataGen(setting = setting, n = n, p = p, beta = beta, rho = rho)
   	X_sp = data$X
    Y_sp = data$Y
    
    fit.proj <- lasso.proj(X_sp, Y_sp-betaNull*X_sp[,which.covariate], standardize = TRUE, parallel = TRUE, ncores = 40)
      
    pval[s,] = fit.proj$pval
       
    if(s %% 100 == 0){  cat("Now computing:", s, "\n")  }
    
  }

  count = 1
  for(j in which.covariate){

    proj.power[count,1] = mean(pval[,j] < 0.2)
    proj.power[count,2] = mean(pval[,j] < 0.1)
    proj.power[count,3] = mean(pval[,j] < 0.05)
    proj.power[count,4] = mean(pval[,j] < 0.01)
  

    count = count + 1
  }
  
  
  return(proj.power)

}  

##############################################################


################ t-test #########################
TTest.Power = function(n = 100, p = 1000, beta, rho, iter = 500, setting = 'dep', which.covariate, betaNull){
  #Return the power of de-sparsified under different settings    
  # Args:
  # setting: different settings, check 'pathwise_simu_setting.R for details'  
  #   rho: related to dependent design setting
  # n,p,beta : sample size, features, coefficients
  #   iter : # of iterations 
  #
  # Return:
  # A matrix of len(which.covariate) * 4 : Simulated power for j-th covariate
  #  
  
  if(length(betaNull) > 1){stop("now only support compute power of 1 coefficients")}
  
  proj.power = matrix(0,length(which.covariate),4)
  
  pval = matrix(NA, iter, p)
  
  for(s in 1:iter){
    
    
    data = dataGen(setting = setting, n = n, p = p, beta = beta, rho = rho)
    
    
    
    X_sp = data$X
    Y_sp = data$Y
    
    
    pval[s,] = (summary(lm(y~., 
                           data = data.frame(y = Y_sp-betaNull*X_sp[,which.covariate], x = X_sp)))$coefficients)[which.covariate+1, 4] 
    
    #fit.proj <- lasso.proj(X_sp, Y_sp-betaNull*X_sp[,which.covariate], standardize = TRUE, parallel = TRUE, ncores = 40)
    #pval[s,] = fit.proj$pval
    
    
    if(s %% 100 == 0){  cat("Now computing:", s, "\n")  }
    
  }
  
  count = 1
  for(j in which.covariate){
    
    proj.power[count,1] = mean(pval[,j] < 0.2)
    proj.power[count,2] = mean(pval[,j] < 0.1)
    proj.power[count,3] = mean(pval[,j] < 0.05)
    proj.power[count,4] = mean(pval[,j] < 0.01)
    
    
    count = count + 1
  }
  
  
  return(proj.power)
  
}  
#################################################

################# F-test #########################

FTest.Power = function(n = 100, p = 1000, beta, rho, iter = 500, setting = 'dep', which.covariate = 1, betaNull = 1){
  #Return the power of de-sparsified under different settings    
  # Args:
  # setting: different settings, check 'pathwise_simu_setting.R for details'  
  #   rho: related to dependent design setting
  # n,p,beta : sample size, features, coefficients
  #   iter : # of iterations 
  #
  # Return:
  # A matrix of len(which.covariate) * 4 : Simulated power for j-th covariate
  #  
  
  #### specified here
  betaNull = 1; which.covariate=1;
  ####
  
  if(length(betaNull) > 1){stop("now only support compute power of 1 coefficients")}
  
  proj.power = matrix(0,length(which.covariate),4)
  
  pval = matrix(NA, iter, p)
  
  for(s in 1:iter){
    
    
    data = dataGen(setting = setting, n = n, p = p, beta = beta, rho = rho)
    
    
    
    X_sp = data$X
    Y_sp = data$Y
    
    
    lm_model = lm(y~., data = data.frame(y=Y_sp-betaNull*X_sp[,which.covariate], x=X_sp))
    
    pval[s,] = linearHypothesis(lm_model, c('x.1=0', 'x.11=0', 'x.12=0'))[[6]][2]
    
    
    #pval[s,] = (summary(lm(y~., 
    #                       data = data.frame(y = Y_sp-betaNull*X_sp[,which.covariate], x = X_sp)))$coefficients)[which.covariate+1, 4] 
    #fit.proj <- lasso.proj(X_sp, Y_sp-betaNull*X_sp[,which.covariate], standardize = TRUE, parallel = TRUE, ncores = 40)
    #pval[s,] = fit.proj$pval
    
    
    if(s %% 100 == 0){  cat("Now computing:", s, "\n")  }
    
  }
  
  count = 1
  for(j in which.covariate){
    
    proj.power[count,1] = mean(pval[,j] < 0.2)
    proj.power[count,2] = mean(pval[,j] < 0.1)
    proj.power[count,3] = mean(pval[,j] < 0.05)
    proj.power[count,4] = mean(pval[,j] < 0.01)
    
    
    count = count + 1
  }
  
  
  return(proj.power)
  
}  

##########################################








