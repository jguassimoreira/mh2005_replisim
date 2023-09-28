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
  
  if (param_2_summarize == "sigma2") { #start with sigma2
    
    sims = sapply(simList, function(x) x[[2]][4,4]) #extract the simulated values from the simList (note we're extracting the variances)
    
    bias = 100 * sims/1; pct_bias = bias - 100
    mean_pct_bias = mean(pct_bias) #compute mean percent bias
    
    hist = ggplot(as.data.frame(pct_bias), aes(x=pct_bias)) +
      geom_histogram(aes(y=..density..), binwidth=1) + geom_vline(aes(xintercept=mean(pct_bias)), color = 'black', linetype="dashed", size = 1.25) +
      geom_density(alpha = 0.2, fill = "green") #create histogram
    
    ci_mat = sapply(simList, function(x) x[[3]][4,]) #extract matrix with CIs for sigma2
    
    coverage = sapply(1:1000, function(x) dplyr::between(1, ci_mat[1,x]^2, ci_mat[2,x]^2)) #check whether each sim is covered
    
    pct_noncoverage = 1 - mean(coverage) #calculate percent noncoverage
    
  } else if (param_2_summarize == 'intercept') {
    
    sims = sapply(simList, function(x) x[[1]][1,1]) #extract the simulated values from the simList
    
    bias = 100 * sims/1; pct_bias = bias - 100
    mean_pct_bias = mean(pct_bias) #compute mean percent bias
    
    hist = ggplot(as.data.frame(pct_bias), aes(x=pct_bias)) +
      geom_histogram(aes(y=..density..), binwidth=1) + geom_vline(aes(xintercept=mean(pct_bias)), color = 'black', linetype="dashed", size = 1.25) +
      geom_density(alpha = 0.2, fill = "green") #create histogram
    
    
    se = sapply(simList, function(x) x[[1]][1,2]) #extract SE for the intercept
    
    ci_mat = matrix(data = c(sims - (2*se), sims + (2*se)), nrow = 2, ncol = 1000, byrow = T) #create CI matrix based on SE and sim values
    
    coverage = sapply(1:1000, function(x) dplyr::between(1, ci_mat[1,x], ci_mat[2,x])) #check whether each sim is covered
    
    pct_noncoverage = 1 - mean(coverage) #calculate percent noncoverage
    
  } else if (param_2_summarize == 'b_x') {
    
    sims = sapply(simList, function(x) x[[1]][2,1]) #extract the simulated values from the simList
    
    bias = 100 * sims/.25; pct_bias = bias - 100
    mean_pct_bias = mean(pct_bias) #compute mean percent bias
    
    hist = ggplot(as.data.frame(pct_bias), aes(x=pct_bias)) +
      geom_histogram(aes(y=..density..), binwidth=1) + geom_vline(aes(xintercept=mean(pct_bias)), color = 'black', linetype="dashed", size = 1.25) +
      geom_density(alpha = 0.2, fill = "green") #create histogram
    
    
    se = sapply(simList, function(x) x[[1]][2,2]) #extract SE for the intercept
    
    ci_mat = matrix(data = c(sims - (2*se), sims + (2*se)), nrow = 2, ncol = 1000, byrow = T) #create CI matrix based on SE and sim values
    
    coverage = sapply(1:1000, function(x) dplyr::between(0.3, ci_mat[1,x], ci_mat[2,x])) #check whether each sim is covered
    
    pct_noncoverage = 1 - mean(coverage) #calculate percent noncoverage
    
  } else if (param_2_summarize == 'b_z') {
    
    sims = sapply(simList, function(x) x[[1]][3,1]) #extract the simulated values from the simList
    
    bias = 100 * sims/.25; pct_bias = bias - 100
    mean_pct_bias = mean(pct_bias) #compute mean percent bias
    
    hist = ggplot(as.data.frame(pct_bias), aes(x=pct_bias)) +
      geom_histogram(aes(y=..density..), binwidth=1) + geom_vline(aes(xintercept=mean(pct_bias)), color = 'black', linetype="dashed", size = 1.25) +
      geom_density(alpha = 0.2, fill = "green") #create histogram
    
    
    se = sapply(simList, function(x) x[[1]][3,2]) #extract SE for the intercept
    
    ci_mat = matrix(data = c(sims - (2*se), sims + (2*se)), nrow = 2, ncol = 1000, byrow = T) #create CI matrix based on SE and sim values
    
    coverage = sapply(1:1000, function(x) dplyr::between(0.3, ci_mat[1,x], ci_mat[2,x])) #check whether each sim is covered
    
    pct_noncoverage = 1 - mean(coverage) #calculate percent noncoverage
    
  } else if (param_2_summarize == 'b_xz') {
    
    sims = sapply(simList, function(x) x[[1]][4,1]) #extract the simulated values from the simList
    
    bias = 100 * sims/.25; pct_bias = bias - 100
    mean_pct_bias = mean(pct_bias) #compute mean percent bias
    
    hist = ggplot(as.data.frame(pct_bias), aes(x=pct_bias)) +
      geom_histogram(aes(y=..density..), binwidth=1) + geom_vline(aes(xintercept=mean(pct_bias)), color = 'black', linetype="dashed", size = 1.25) +
      geom_density(alpha = 0.2, fill = "green") #create histogram
    
    
    se = sapply(simList, function(x) x[[1]][4,2]) #extract SE for the intercept
    
    ci_mat = matrix(data = c(sims - (2*se), sims + (2*se)), nrow = 2, ncol = 1000, byrow = T) #create CI matrix based on SE and sim values
    
    coverage = sapply(1:1000, function(x) dplyr::between(0.3, ci_mat[1,x], ci_mat[2,x])) #check whether each sim is covered
    
    pct_noncoverage = 1 - mean(coverage) #calculate percent noncoverage
    
    
  } else if (param_2_summarize == 'int_var') {
    
    sims = sapply(simList, function(x) x[[2]][1,4]) #extract the simulated values from the simList
    
    bias = 100 * sims/calc_var_comp_legacy(icc, 1)[1]; pct_bias = bias - 100
    mean_pct_bias = mean(pct_bias) #compute mean percent bias
    
    hist = ggplot(as.data.frame(pct_bias), aes(x=pct_bias)) +
      geom_histogram(aes(y=..density..), binwidth=1) + geom_vline(aes(xintercept=mean(pct_bias)), color = 'black', linetype="dashed", size = 1.25) +
      geom_density(alpha = 0.2, fill = "green") #create histogram
    
    
    ci_mat = sapply(simList, function(x) x[[3]][1,]) #extract CIs 
    
    coverage = sapply(1:1000, function(x) dplyr::between(calc_var_comp_legacy(icc, 1)[1], ci_mat[1,x]^2, ci_mat[2,x]^2)) #check whether each sim is covered
    
    pct_noncoverage = 1 - mean(coverage) #calculate percent noncoverage
    
    
  } else if (param_2_summarize == 'x_var') {
    
    sims = sapply(simList, function(x) x[[2]][2,4]) #extract the simulated values from the simList
    
    bias = 100 * sims/calc_var_comp_legacy(icc, 1)[1]; pct_bias = bias - 100
    mean_pct_bias = mean(pct_bias) #compute mean percent bias
    
    hist = ggplot(as.data.frame(pct_bias), aes(x=pct_bias)) +
      geom_histogram(aes(y=..density..), binwidth=1) + geom_vline(aes(xintercept=mean(pct_bias)), color = 'black', linetype="dashed", size = 1.25) +
      geom_density(alpha = 0.2, fill = "green") #create histogram
    
    
    ci_mat = sapply(simList, function(x) x[[3]][3,]) #extract CIs 
    
    coverage = sapply(1:1000, function(x) dplyr::between(calc_var_comp_legacy(icc, 1)[1], ci_mat[1,x]^2, ci_mat[2,x]^2)) #check whether each sim is covered
    
    pct_noncoverage = 1 - mean(coverage) #calculate percent noncoverage
    
  }
  
  outList = list(mean_pct_bias, hist, coverage, pct_noncoverage)
  names(outList) = c('mean_pct_bias', 'hist', 'coverage', 'pct_noncoverage')
  
  return(outList) #return the tau matrix
  
}
  