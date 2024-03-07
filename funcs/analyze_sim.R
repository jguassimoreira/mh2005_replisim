##############################################
## function to analyze output of simulation ##

#intputs: param_2_analyze - parameter you would like analyzed
#         conds - vector of strings indicating the conditions from the simulation
#         simList - list of simulation output, broken down by condition

#outputs: TBD FILL THIS IN, JOAO

analyze_sim = function(simList, param_2_summarize, conds) {
  
  ## First we're going to start with coverage
  
  #extract icc conditions
  icc_conds = lapply(conds, function(x) parse_conds(x)[[1]])
  
  #extract coverage outcomes for each replication across all conditions
  sim_datasets = lapply(1:length(conds), function(x) cbind(summarize_sim(simList[[x]], param_2_summarize, icc_conds[[x]])[[3]],
                                                                   icc_conds[[x]], 
                                                                   parse_conds(conds[x])[[2]],
                                                                   parse_conds(conds[x])[[3]]))
  
  #since output of lapply is a list, bind these together into 1 big dataframe
  sim_datasets = as.data.frame(do.call(rbind, sim_datasets))
  
  #name the columns in the dataset
  names(sim_datasets) = c("coverage", "icc", "ngrps", "grpsize")
  
  #make coverage values numeric, flip direction (so 1 means *non*coverage)
  sim_datasets$coverage = ifelse(as.numeric(as.factor(sim_datasets$coverage)) == 2, 0, 1)
  #turn the rest of the features/variables into factors
  sim_datasets$icc = as.factor(sim_datasets$icc)
  sim_datasets$ngrps = as.factor(sim_datasets$ngrps)
  sim_datasets$grpsize = as.factor(sim_datasets$grpsize)
  
  #logit models
  logitmodel_all = glm(coverage ~ icc + ngrps + grpsize, data = sim_datasets, family = "binomial") #all factors entered simultaneously
  
  logitmodel_icc = glm(coverage ~ icc, data = sim_datasets, family = "binomial") #just the ICC factor entered
  
  logitmodel_ngrps = glm(coverage ~ ngrps, data = sim_datasets, family = "binomial") #just the ngrps factor entered
  
  logitmodel_grpsize = glm(coverage ~ grpsize, data = sim_datasets, family = "binomial") #just the grpsize factor entered
  
  logitmodel_interaction = glm(coverage ~ icc*ngrps + icc*grpsize + ngrps*grpsize, data = sim_datasets, family = "binomial") #all two interactions entered
  
  #wald tests of the three factors
  
  #first the models with only 1 factor entered
  #waldicc = wald.test(b = coef(logitmodel_icc), Sigma = vcov(logitmodel_icc), Terms = 2:3)
  #waldngrps = wald.test(b = coef(logitmodel_ngrps), Sigma = vcov(logitmodel_ngrps), Terms = 2:3)
  #waldgrpsize = wald.test(b = coef(logitmodel_grpsize), Sigma = vcov(logitmodel_grpsize), Terms = 2:3)
  waldicc = anova(logitmodel_icc, test = "Chisq")
  waldngrps = anova(logitmodel_ngrps, test = "Chisq")
  waldgrpsize = anova(logitmodel_grpsize, test = "Chisq")
  
  #next the models with all the factors entered simultaneously
  waldicc_all = wald.test(b = coef(logitmodel_all), Sigma = vcov(logitmodel_all), Terms = 2:3)
  waldngrps_all = wald.test(b = coef(logitmodel_all), Sigma = vcov(logitmodel_all), Terms = 4:5)
  waldgrpsize_all = wald.test(b = coef(logitmodel_all), Sigma = vcov(logitmodel_all), Terms = 6:7)
  
  #now do the ones for the interactions
  waldicc_ngrps = wald.test(b = coef(logitmodel_interaction), Sigma = vcov(logitmodel_interaction), Terms = 8:11)
  waldicc_grpsize = wald.test(b = coef(logitmodel_interaction), Sigma = vcov(logitmodel_interaction), Terms = 12:15)
  waldngrps_grpsize = wald.test(b = coef(logitmodel_interaction), Sigma = vcov(logitmodel_interaction), Terms = 16:19)
  
  #aggregate mean coverages broken down by the three factors, and then compute aggregate mean coverages for the fully crossed design (i.e., for each cell)
  agg_icc = aggregate(sim_datasets$coverage, by = list(sim_datasets$icc), FUN = mean, na.rm = T)
  agg_ngrps = aggregate(sim_datasets$coverage, by = list(sim_datasets$ngrps), FUN = mean, na.rm = T)
  agg_grpsize = aggregate(sim_datasets$coverage, by = list(sim_datasets$grpsize), FUN = mean, na.rm = T)
  
  agg_crossed = aggregate(sim_datasets$coverage, by = list(sim_datasets$icc, sim_datasets$ngrps, sim_datasets$grpsize), FUN = mean, na.rm = T)
  
  #put outputs in lists, nest into a final list, return
  logitList = list(logitmodel_all, logitmodel_icc, logitmodel_ngrps, logitmodel_grpsize, logitmodel_interaction)
  waldList = list(waldicc, waldngrps, waldgrpsize, waldicc_all, waldngrps_all, waldgrpsize_all, waldicc_ngrps, waldicc_grpsize, waldngrps_grpsize)
  names(waldList) = c("icc", 'ngrps', 'grpsize', "icc_all", 'ngrps_all', 'grpsize_all', 'icc_ngrps', 'icc_grpsize', 'ngrps_grpsize')
  aggList = list(agg_icc, agg_ngrps, agg_grpsize, agg_crossed); names(aggList) = c('icc', 'ngrps', 'grpsize', 'crossed')
  
  ##Now we'll move on to percent bias
  sim_bias = lapply(1:length(conds), function(x) cbind(summarize_sim(simList[[x]], param_2_summarize, icc_conds[[x]])[[1]],
                                                       icc_conds[[x]], 
                                                       parse_conds(conds[x])[[2]],
                                                       parse_conds(conds[x])[[3]]))
  sim_bias = data.frame(do.call('rbind', sim_bias))
  names(sim_bias) = c("bias", "icc", "ngrps", "grpsize")
  sim_bias$bias = as.numeric(sim_bias$bias)
  icc_bias = aggregate(bias ~ icc, data = sim_bias[,c('bias', 'icc')], mean, na.rm = T); names(icc_bias) = c('cond_lev', 'bias')
  ngrps_bias = aggregate(bias ~ ngrps, data = sim_bias[,c('bias', 'ngrps')], mean, na.rm = T); names(ngrps_bias) = c('cond_lev', 'bias')
  grpsize_bias = aggregate(bias ~ grpsize, data = sim_bias[,c('bias', 'grpsize')], mean, na.rm = T); names(grpsize_bias) = c('cond_lev', 'bias')
  
  bias_marg = rbind(icc_bias, ngrps_bias, grpsize_bias)
  
  outList = list(bias_marg, logitList, waldList, aggList); names(outList) = c('bias', 'logitmodel', 'wald_tests', 'agg_coverage')
  
  return(outList)
  
}