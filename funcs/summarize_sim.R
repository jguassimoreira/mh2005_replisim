################################################
## function to summarize output of simulation ##

#intputs: simList - list of simulation output, broken down by condition
#         param_2_summarize - parameter you would like summarized
#         icc - intraclass correlation coefficient used to simulate the dataset you want summarized

#outputs: A list containing the mean percent bias in parameter estimation (scalar) across reps for
#         a given condition, a histogram+density plot of the bias values across replications
#         a vector (coverage) indicating whether 95% CI for a given rep contained true param value (set by exp)
#         TRUE indicates 95% for that replication did indeed contain the true parameter value
#         pct_noncoverage is the proportion of replications whose CIs did not contain true value (scalar)

summarize_sim = function(simList, param_2_summarize, icc) {
  
  true_params = list('sigma2' = 0.5, 'intercept' = 1, 'b_x' = 0.3,
                     'b_z' = 0.3, 'b_xz' = 0.3, 'int_var' = calc_var_comp_legacy(icc, .5)[1], 
                     'x_var' = calc_var_comp_legacy(icc, .5)[1])
  
  if (param_2_summarize == "sigma2") {
    
    #extract the estimated values fit to the model from the simulated dataset
    sims = sapply(simList, function(x) x$varComp[x$varComp$grp == 'Residual', 'vcov'])
    
    #get the confidence interval values
    #confidence intervals for variances are on SD scale, need to square to put back into variance
    ci_mat = sapply(simList, function(x) x$confInt['.sigma',]^2)
    
  } else if (param_2_summarize == 'intercept') {
    
    #extract the estimated values fit to the model from the simulated dataset
    sims = sapply(simList, function(x) x$fx_coefs['(Intercept)', 'Estimate'])
    
    #get the confidence interval values
    ci_mat = sapply(simList, function(x) x$confInt['(Intercept)',])
    
  } else if (param_2_summarize == 'b_x') {
    
    #extract the estimated values fit to the model from the simulated dataset
    sims = sapply(simList, function(x) x$fx_coefs['x', 'Estimate'])
    
    #get the confidence interval values
    ci_mat = sapply(simList, function(x) x$confInt['x',])
    
  } else if (param_2_summarize == 'b_z') {
    
    #extract the estimated values fit to the model from the simulated dataset
    sims = sapply(simList, function(x) x$fx_coefs['z', 'Estimate'])
    
    #get the confidence interval values
    ci_mat = sapply(simList, function(x) x$confInt['z',])
    
  } else if (param_2_summarize == 'b_xz') {
    
    #extract the estimated values fit to the model from the simulated dataset
    sims = sapply(simList, function(x) x$fx_coefs['x:z', 'Estimate'])
    
    #get the confidence interval values
    ci_mat = sapply(simList, function(x) x$confInt['x:z',])
    
  } else if (param_2_summarize == 'int_var') {
    
    #extract the estimated values fit to the model from the simulated dataset
    sims = sapply(simList, function(x) x$varComp[x$varComp$var1 == '(Intercept)' & is.na(x$varComp$var2), 'vcov'][1])
    
    #get the confidence interval values
    #confidence intervals for variances are on SD scale, need to square to put back into variance
    ci_mat = sapply(simList, function(x) x$confInt['.sig01',]^2)
    
  } else if (param_2_summarize == 'x_var') {
    
    #extract the estimated values fit to the model from the simulated dataset
    sims = sapply(simList, function(x) x$varComp[x$varComp$var1 == 'x' & is.na(x$varComp$var2), 'vcov'][1])
    
    #get the confidence interval values
    #confidence intervals for variances are on SD scale, need to square to put back into variance
    ci_mat = sapply(simList, function(x) x$confInt['.sig03',]^2)
    
  }
  
  #compute percent bias
  bias = 100 * sims/true_params[[param_2_summarize]]; pct_bias = bias - 100 
  #get mean bias
  mean_pct_bias = mean(pct_bias) 
  
  #create histogram
  hist = ggplot(as.data.frame(pct_bias), aes(x=pct_bias)) +
    geom_histogram(aes(y=..density..), binwidth=1) + geom_vline(aes(xintercept=mean(pct_bias)), color = 'black', linetype="dashed", size = 1.25) +
    geom_density(alpha = 0.2, fill = "green") 
  
  #get coverage values
  coverage = sapply(1:1000, function(x) dplyr::between(true_params[[param_2_summarize]], ci_mat[1,x], ci_mat[2,x]))
  
  #calculate percent non-coverage
  pct_noncoverage = 1 - mean(coverage, na.rm = T) #calculate percent noncoverage
  
  #get non-covergence
  non_converged = sum(sapply(simList, function(x) x$w))
  
  #create list of outputs
  outList = list(mean_pct_bias, hist, coverage, pct_noncoverage, non_converged)
  names(outList) = c('mean_pct_bias', 'hist', 'coverage', 'pct_noncoverage', 'non_converged')
  
  return(outList) #return the list of summarized quantities
  
}
  