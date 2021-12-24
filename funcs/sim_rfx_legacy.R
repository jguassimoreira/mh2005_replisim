#############################################################################
## function to simulate random effects for a simulated dataset via M&H2005 ##

#intputs: N_l2 - number of level 2 units
#         mu - vector of random effect means, usually set to zero
#         varcovar - variance/covariance matrix for random effects

#outputs: A matrix containing random effects for each level 2 unit (units - rows; rfx - columns)


sim_rfx_legacy = function(N_l2, muVec, varcovar){
  
  #follow Maas & Hox's 2005 specifications
  rfx = MASS::mvrnorm(n = N_l2, #number of L2 units
                      muVec = mu, #vector of rfx means, usually zero
                      Sigma = varcovar) #variance components
  
  return(rfx)
  
}