##########################################################
## function to run MLM simulations that extend M&H 2005 ##

#inputs: N_l2 - level 2 sample size
#         n_l1 - level 1 sample size
#         mod_obj - mlm power model object, used as part of data generation
#         imbalance - ratio of size between largest and smallest groups. Largest group always has full n_l1 size. 

#outputs: Matrix of fixed effect estimates, standard errors, t stats and p-values

##Part of the code in this function was adapted from the mlmpower package (https://github.com/bkeller2/mlmpower)
##See Enders, Keller, & Woller 2023, Psych Methods (doi: https://doi.org/10.1037/met0000614)

sim_mlm_extension = function(N_l2, n_l1, mod_obj, imbalance) { 
  
  ### Generate data
  
  n_within = n_l1 #set within sample size
  n_between = N_l2 #set between sample size
  
  p = make_parameters(mod_obj) #get model parameters from mlm_power model object
  
  n_within_imbalanced = ceiling(runif(n_between, n_within/imbalance, n_within)) #sample imbalanced group sizes from uniform dist
  
  N = sum(n_within_imbalanced) #Total sample size (was: n_within * n_between)
  l1 = length(p$mean_X) #number of level 1 predictors
  l2 = length(p$mean_Z) #number of level 2 predictors
  
  # Generate ID variable
  `_id` = matrix(rep(1:n_between, time = n_within_imbalanced), ncol = 1) #was: seq_len(n_between) %x% matrix(1, n_within)
  
  # Generate within and between variables (for X matrix)
  X_w = rmvnorm_nomean(N, p$phi_w) #x - within
  X_b = rmvnorm(N, c(p$mean_X, rep(0, l2)), p$phi_b)[`_id`, , drop = F] #x - between (includes mean of within of var and L2 var)
  
  # Create all interactions
  inter = do.call('cbind', apply(
    X_w, 2,
    \(.). * X_b[, seq_len(l2) + l1, drop = F],
    simplify = F
  ))
  
  # Create predictors matrix
  X = cbind(1, X_w, inter, X_b)
  
  # Generate level-1 residuals
  e_ij <- rnorm(N, 0, sqrt(p$var_e))
  
  # Generate level-2 residuals
  u_j <- rmvnorm_nomean(n_between, p$tau)[`_id`, , drop = F]
  
  # Generate Outcome
  Y <- X %*% p$gammas + rowSums(X[,seq_len(l1 + 1), drop = F] * u_j) + e_ij
  
  # Compile into dataframe
  d <- data.frame(
    `_id`,
    Y,
    X_w + X_b[ , seq_len(l1), drop = F],
    (X_b[ , seq_len(l2) + l1, drop = F]
     # Add means for Z
     + matrix(
       p$mean_Z,
       nrow = N,
       ncol = l2,
       byrow = T
     )
    )
  )
  
  #add column names
  names(d) <- c(
    '_id',
    names(p$mean_Y),
    names(which(vapply(mh$predictors, levels, numeric(1L)) == 1)),
    names(which(vapply(mh$predictors, levels, numeric(1L)) == 2))
  )
  
  ### Fit the model
  
  #set the model string
  
  test = 'y ~ x1 + z1 + x1*z1 + (1+x1|`_id`)'
  
  mod_string = sprintf("y ~ %s + %s",
                       paste0(names(mod_obj$l1), "*z1", collapse = " + "),
                       paste0("1+", names(mod_obj$l1), collapse = "+"))
  
  mod = lmerTest::lmer(y ~ x + z + x*z + (1+x|subID), data = dat, REML = T, lme4::lmerControl(check.conv.singular = lme4::.makeCC(action = "ignore",  tol = 1e-4)))
  
  lmerTest::lmer(test, data = d, REML = T, lme4::lmerControl(check.conv.singular = lme4::.makeCC(action = "ignore",  tol = 1e-4)))
  }