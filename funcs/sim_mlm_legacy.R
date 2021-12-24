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
  
  dat = agument_dat_legacy(n_l1, N_l2, xdat, zdat) #augment predictor design matrix (i.e., repeat observations if sample size combo isn't the smallest)
  
  dat$xz = dat$x * dat$z # compute cross product for interaction term 
  
  #simulate the 'true' y (or mean/expected Y given many iterations of the same experiment)
  ymn = 1 + (0.3*dat$z) + (0.3*dat$x) + (0.3*dat$xz) + rep(rfx[,1], each = n_l1) + (rep(rfx[,2], each = n_l1)*dat$x) 
  #according to M&H, intercept (gamma00) should be equal to 1, and all other coefficients equal to 0.3 (gamma01,gamma10,gamma11)
  #first term is the intercept (1, gamma00), next three terms are fixed effect slopes (gamma01, gamma10, gamma11)
  #the following term (2nd to last) is random intercept (technically residual portion of random intercept)
  #last term is interaction between x and random slope (residual portion of random slope)
  
  l1errors = rnorm(n_l1*N_l2, 0 , sqrt(sigma2)) #simulate the l1 error
  #technically we don't need to take square of sigma2 since M&H set to 1, but good to code to allow for future flexibility
  
  dat$y = ymn + l1errors #add y to the data frame as sum of expected Y given predictors&random fx + l1 error
  
  mod = lmer(y ~ x + z + x*z + (1+x|subID), data = dat, REML = T) #run the multilevel model with REML (as specified by M&H)
  
  #now it's time to extract the stuff we need and return it in a list
  
  vc = lme4::VarCorr(mod)#extract estimated tau matrix from model
  varcomp = c(vc[[1]][1], vc[[1]][2], vc[[1]][3], vc[[1]][4])
  
  fx = coef(summary(mod)) #extract estimated fixed effects
  
  modsigma2 = (mod@sigma)^2 #extract estimated sigma, square it
  
  #the last we'll pull out is the icc of the x variable, to better understand consequences of M&H's predictor generation scheme 
  xiccmod = lmer(x ~ 1 + (1|subID), data = dat, REML = T) #run empty model with x variable as the DV
  xicc = lme4::VarCorr(xiccmod)[[1]][1] / (lme4::VarCorr(xiccmod)[[1]][1] + xiccmod@sigma^2) #pull out ICC
  
  return(list(fx, modsigma2, varcomp, xicc)) #return these items as a list
  
}

