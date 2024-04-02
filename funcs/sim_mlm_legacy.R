###########################################################
## function to run MLM simulations according to M&H 2005 ##

#inputs: N_l2 - level 2 sample size
#         n_l1 - level 1 sample size
#         icc - intraclass coefficient
#         sigma2 - level 1 error variance (analagous to MSE in OLS regression)
#         xdat - fixed level 1 (within cluster) data
#         zdat - fixed level 2 (between cluster) data

#outputs: A list containing the following elements
#         Matrix of fixed effect estimates, standard errors, t stats and p-values
#         Scalar value denoting sigma2

sim_mlm_legacy = function(N_l2, n_l1, icc, sigma2, xdat, zdat) {

  run = TRUE
  
  while (run == TRUE) {
    
    w = 0
    
    varComp = calc_var_comp_legacy(icc, sigma2) #calculate variance components (tau00 and tau11)
    
    rfx = sim_rfx_legacy(N_l2, c(0,0), varComp) #use variance components to simulate random effects
    
    dat = augment_dat_legacy(n_l1, N_l2, xdat, zdat) #augment predictor design matrix (i.e., repeat observations if sample size combo isn't the smallest)
    
    dat$xz = dat$x * dat$z # compute cross product for interaction term 
    
    #simulate the 'true' y (or mean/expected Y given many iterations of the same experiment)
    ymn = 1 + (0.3*dat$z) + (0.3*dat$x) + (0.3*dat$xz) + rep(rfx[,1], each = n_l1) + (rep(rfx[,2], each = n_l1)*dat$x) 
    #according to M&H, intercept (gamma00) should be equal to 1, and all other coefficients equal to 0.3 (gamma01,gamma10,gamma11)
    #first term is the intercept (1, gamma00), next three terms are fixed effect slopes (gamma01, gamma10, gamma11)
    #the following term (2nd to last) is random intercept (technically residual portion of random intercept)
    #last term is interaction between x and random slope (residual portion of random slope)
    
    l1errors = rnorm(n_l1*N_l2, 0 , sqrt(sigma2)) #simulate the l1 error
    
    dat$y = ymn + l1errors #add y to the data frame as sum of expected Y given predictors&random fx + l1 error
    
    mod = lmerTest::lmer(y ~ x + z + x*z + (1+x|subID), data = dat, REML = T, lme4::lmerControl(check.conv.singular = lme4::.makeCC(action = "ignore",  tol = 1e-4))) #run the multilevel model with REML (as specified by M&H)
    
    #do some checks
    warn = is_warning_generated(mod) #check to make sure model converged
    sing = isSingular(mod) #check if the model has singular fit
    
    #if neither error is raised, we're good to go!
    if (warn == FALSE & sing == FALSE) {
      
      run = FALSE
      
    } else {
      
      #if either error is raised, we keep simulating data and use w to count the number of extra sims neededs
      run = TRUE
      w = w+1
      
    }
  }
    
  #once we've made it here, we assume that the model converged without error
  #now it's time to extract the stuff we need and return it in a list
    
  varcomp = data.frame(lme4::VarCorr(mod)) #save the variance components
  #conf_ints = lme4::confint.merMod(mod, method ="profile", oldNames = F, quiet = T) #get the confidence intervals for all the model parameters
  conf_ints_var = compute_wald_ci_varcomp(mod) #compute confidence intervals for variance components
    
  conf_ints_oth = confint(mod, method = "Wald", quiet = T) #compute confidence intervals for all other params
    
  conf_ints = rbind(conf_ints_oth[c('(Intercept)', 'x', 'z', 'x:z'),], conf_ints_var) #combine the confidence intervals we need
  #taking rows 1:3 of conf_ints_var gets just the CIs of the tau matrix elements
    
  fx = coef(summary(mod)) #extract estimated fixed effects
  
  return(list(fx, varcomp, conf_ints, w)) #return these items as a list
  
}

