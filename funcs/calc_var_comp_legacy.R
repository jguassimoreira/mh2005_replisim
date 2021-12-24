###################################################################################
## function to calculate variance components for a simulated dataset via M&H2005 ##

#intputs: icc - intraclass correlation coefficient (for y variable)
#         sigma2 - level 1 error variance, i.e., average squared within cluster error (MSE analog from OLS)

#outputs: The tau matrix containing the variance components (under constraints specified by M&H)


calc_var_comp_legacy = function(icc, sigma2){
  
  tau00 = (icc*sigma2) / (1-icc) #define residual random intercept variance; derived from eqn 6 in M&H
  tau11 = tau00 #define residual random slope variance by setting equal to tau00; justified by M&H on page 89, 1st para
  tau01 = 0 #setting covariance equal to zero, as specified by M&H
  
  tauMat = matrix(c(tau00, tau01, tau01, tau11), nrow = 2, byrow = F) #put variance components into a matrix (tau matrix)
  
  return(tauMat) #return the tau matrix
  
}
  