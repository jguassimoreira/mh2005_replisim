###########################################################
## function to run MLM simulations according to M&H 2005 ##

#intputs: N_l2 - level 2 sample size
#         n_l1 - level 1 sample size
#         icc - intraclass coefficient
#         sigma2 - level 1 error variance (analagous to MSE in OLS regression)
#         xdat - fixed level 1 (within cluster) data
#         zdat - fixed level 2 (between cluster) data

#outputs: A list containing the following elements
#         Matrix of fixed effect estimates, standard errors, t stats and p-values
#         Scalar value denoting sigma2


sim_mlm_legacy = function(N_l2, n_l1, icc, sigma2, xdat, zdat) {

  varComp = calc_var_comp_legacy(icc, sigma2) #calculate variance components (tau00 and tau11)
  
  rfx = sim_rfx_legacy(N_l2, c(0,0), varComp) #use variance components to simulate random effects
  
  
  
}

