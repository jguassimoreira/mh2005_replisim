##############################################
## function to analyze output of simulation ##

#intputs: param_2_analyze - parameter you would like analyzed
#         conds - vector of strings indicating the conditions from the simulation
#         simList - list of simulation output, broken down by condition

#outputs: TBD FILL THIS IN, JOAO

analyze_sim = function(simList, param_2_summarize, conds) {
  
  
  icc_conds = lapply(1:length(conds), function(x) ifelse(strsplit(strsplit(conds[x], '/')[[1]][7], "_")[[1]][1] == 'icc1', 0.1,
                                                         ifelse(strsplit(strsplit(conds[x], '/')[[1]][7], "_")[[1]][1] == 'icc2', 0.2, 0.3)))
  
  #extract coverage outcomes for each replication across all conditions
  sim_datasets = lapply(1:length(conds), function(x) cbind(summarize_sim(simList[[x]], param_2_summarize, icc_conds[[x]])[[3]], 
                                                                   rep(strsplit(strsplit(conds[x], '/')[[1]][7], "_")[[1]][1], 1000), #yep
                                                                   rep(strsplit(strsplit(conds[x], '/')[[1]][7], "_")[[1]][2], 1000),
                                                                   rep(strsplit(strsplit(conds[x], '/')[[1]][7], "_")[[1]][3], 1000)))
  
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
  waldicc = wald.test(b = coef(logitmodel_icc), Sigma = vcov(logitmodel_icc), Terms = 2)
  waldngrps = wald.test(b = coef(logitmodel_ngrps), Sigma = vcov(logitmodel_ngrps), Terms = 2)
  waldgrpsize = wald.test(b = coef(logitmodel_grpsize), Sigma = vcov(logitmodel_grpsize), Terms = 2)
  
  #next the models with all the factors entered simultaneously
  waldicc_all = wald.test(b = coef(logitmodel_all), Sigma = vcov(logitmodel_all), Terms = 2:3)
  waldngrps_all = wald.test(b = coef(logitmodel_all), Sigma = vcov(logitmodel_all), Terms = 4:5)
  waldgrpsize_all = wald.test(b = coef(logitmodel_all), Sigma = vcov(logitmodel_all), Terms = 6:7)
  
  #now do the ones for the interactions
  waldicc_ngrps = wald.test(b = coef(logitmodel_interaction), Sigma = vcov(logitmodel_interaction), Terms = 8:11)
  waldicc_grpsize = wald.test(b = coef(logitmodel_interaction), Sigma = vcov(logitmodel_interaction), Terms = 12:15)
  waldngrps_grpsize = wald.test(b = coef(logitmodel_interaction), Sigma = vcov(logitmodel_interaction), Terms = 16:19)
  
  #aggregate mean coverages broken down by the three factors, and then compute aggregate mean coverages for the fully crossed design (i.e., for each cell)
  agg_icc = aggregate(sim_datasets$coverage, by = list(sim_datasets$icc), FUN = mean)
  agg_ngrps = aggregate(sim_datasets$coverage, by = list(sim_datasets$ngrps), FUN = mean)
  agg_grpsize = aggregate(sim_datasets$coverage, by = list(sim_datasets$grpsize), FUN = mean)
  
  agg_crossed = aggregate(sim_datasets$coverage, by = list(sim_datasets$icc, sim_datasets$ngrps, sim_datasets$grpsize), FUN = mean)
  
  #put outputs in lists, nest into a final list, return
  logitList = list(logitmodel_all, logitmodel_icc, logitmodel_ngrps, logitmodel_grpsize, logitmodel_interaction)
  waldList = list(waldicc, waldngrps, waldgrpsize, waldicc_all, waldngrps_all, waldgrpsize_all, waldicc_ngrps, waldicc_grpsize, waldngrps_grpsize)
  names(waldList) = c("icc", 'ngrps', 'grpsize', "icc_all", 'ngrps_all', 'grpsize_all', 'icc_ngrps', 'icc_grpsize', 'ngrps_grpsize')
  aggList = list(agg_icc, agg_ngrps, agg_grpsize, agg_crossed); names(aggList) = c('icc', 'ngrps', 'grpsize', 'crossed')
  
  outList = list(logitList, waldList, aggList); names(outList) = c('logitmodel', 'wald_tests', 'agg_coverage')
  
  return(outList)
  
}