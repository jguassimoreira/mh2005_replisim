#################################################################
## function to run interaction analysis of the replication sim ##

replication_interaction_analysis = function(simList, param_2_analyze, conds) {
  
  outList = list("param" = param_2_analyze) #initalize output list here, we'll append to it as we go along
  
  icc_conds = lapply(conds, function(x) parse_conds(x)[[1]])
  
  #extract coverage outcomes for each replication across all conditions
  sim_datasets = lapply(1:length(conds), function(x) cbind(summarize_sim(simList[[x]], param_2_analyze, icc_conds[[x]])[[3]],
                                                           icc_conds[[x]], 
                                                           parse_conds(conds[x])[[2]],
                                                           parse_conds(conds[x])[[3]]))
  
  #since output of lapply is a list, bind these together into 1 big dataframe
  sim_datasets = as.data.frame(do.call(rbind, sim_datasets))
  
  #name the columns in the dataset
  names(sim_datasets) = c("coverage", "icc", "ngrps", "grpsize")
  
  #create *non* coverage variable as numeric
  sim_datasets$non_coverage = dplyr::case_when(
    sim_datasets$coverage == "TRUE" ~ 0, #true coverage means CIs covered the ground truth param, so we need to flip here
    sim_datasets$coverage == "FALSE" ~ 1 #same as above
  )
  
  #turn the rest of the features/variables into factors
  sim_datasets$icc = as.factor(sim_datasets$icc)
  sim_datasets$ngrps = as.factor(sim_datasets$ngrps)
  sim_datasets$grpsize = as.factor(sim_datasets$grpsize)
  
  #highest level interaction (here 3 way)
  all_int_mod = glm(non_coverage ~ icc*ngrps*grpsize, data = sim_datasets, family = "binomial") 
  
  #compute LRT for highest level model
  all_int_lrt = anova(all_int_mod, test = "LRT")
  
  #get p-value for the highest interaction
  p_val_all_int = all_int_lrt['icc:ngrps:grpsize','Pr(>Chi)']
  
  if (p_val_all_int >= .01) {
    
    #save this message if it's not significant
    all_int_test_msg = 'Not significant'
    outList[["all_int_result"]] = all_int_test_msg
    
  } else {
    
    #if it is significant, spit out the p-value and then test all possible two interactions at each level of the third variable
    #then get the p-values for those guys and save them
    #we'll spit them out and the user can proceed from there as needed
    all_int_test_msg = sprintf('significant, p = %s', p_val_all_int)
    outList[["all_int_result"]] = all_int_test_msg
    
    #test ngrps * grpsize at all possible levels of icc
    ngrps_grpsize_icc1 = glm(non_coverage ~ ngrps*grpsize, data = sim_datasets[sim_datasets$icc == 0.1,], family = "binomial")
    ngrps_grpsize_icc2 = glm(non_coverage ~ ngrps*grpsize, data = sim_datasets[sim_datasets$icc == 0.2,], family = "binomial")
    ngrps_grpsize_icc3 = glm(non_coverage ~ ngrps*grpsize, data = sim_datasets[sim_datasets$icc == 0.3,], family = "binomial")
    
    #test ngrps * icc at all possible levels of grpsize
    ngrps_icc_grpsize5 = glm(non_coverage ~ ngrps*icc, data = sim_datasets[sim_datasets$grpsize == "grpsize5",], family = "binomial")
    ngrps_icc_grpsize30 = glm(non_coverage ~ ngrps*icc, data = sim_datasets[sim_datasets$grpsize == "grpsize30",], family = "binomial")
    ngrps_icc_grpsize50 = glm(non_coverage ~ ngrps*icc, data = sim_datasets[sim_datasets$grpsize == "grpsize50",], family = "binomial")
    
    #test grpsize*icc at all possible levels of ngrps
    grpsize_icc_ngrps30  = glm(non_coverage ~ grpsize*icc, data = sim_datasets[sim_datasets$ngrps == "ngrps30",], family = "binomial")
    grpsize_icc_ngrps50  = glm(non_coverage ~ grpsize*icc, data = sim_datasets[sim_datasets$ngrps == "ngrps50",], family = "binomial")
    grpsize_icc_ngrps100  = glm(non_coverage ~ grpsize*icc, data = sim_datasets[sim_datasets$ngrps == "ngrps100",], family = "binomial")
    
    #compute likelihood ratio tests
    ngrps_grpsize_icc1_lrt = anova(ngrps_grpsize_icc1, test = "LRT")
    ngrps_grpsize_icc2_lrt = anova(ngrps_grpsize_icc2, test = "LRT")
    ngrps_grpsize_icc3_lrt = anova(ngrps_grpsize_icc3, test = "LRT")
    
    ngrps_icc_grpsize5_lrt = anova(ngrps_icc_grpsize5, test = "LRT")
    ngrps_icc_grpsize30_lrt = anova(ngrps_icc_grpsize30, test = "LRT")
    ngrps_icc_grpsize50_lrt = anova(ngrps_icc_grpsize50, test = "LRT")
    
    grpsize_icc_ngrps30_lrt  = anova(grpsize_icc_ngrps30, test = "LRT")
    grpsize_icc_ngrps50_lrt  = anova(grpsize_icc_ngrps50, test = "LRT")
    grpsize_icc_ngrps100_lrt = anova(grpsize_icc_ngrps100, test = "LRT")
    
    #get pvalues
    ngrps_grpsize_icc1_p = ngrps_grpsize_icc1_lrt[,'Pr(>Chi)'][length(ngrps_grpsize_icc1_lrt[,'Pr(>Chi)'])]
    ngrps_grpsize_icc2_p = ngrps_grpsize_icc2_lrt[,'Pr(>Chi)'][length(ngrps_grpsize_icc2_lrt[,'Pr(>Chi)'])]
    ngrps_grpsize_icc3_p = ngrps_grpsize_icc3_lrt[,'Pr(>Chi)'][length(ngrps_grpsize_icc3_lrt[,'Pr(>Chi)'])]
    
    ngrps_icc_grpsize5_p = ngrps_icc_grpsize5_lrt[,'Pr(>Chi)'][length(ngrps_icc_grpsize5_lrt[,'Pr(>Chi)'])]
    ngrps_icc_grpsize30_p = ngrps_icc_grpsize30_lrt[,'Pr(>Chi)'][length(ngrps_icc_grpsize30_lrt[,'Pr(>Chi)'])]
    ngrps_icc_grpsize50_p = ngrps_icc_grpsize50_lrt[,'Pr(>Chi)'][length(ngrps_icc_grpsize50_lrt[,'Pr(>Chi)'])]
    
    grpsize_icc_ngrps30_p = grpsize_icc_ngrps30_lrt[,'Pr(>Chi)'][length(grpsize_icc_ngrps30_lrt[,'Pr(>Chi)'])]
    grpsize_icc_ngrps50_p = grpsize_icc_ngrps50_lrt[,'Pr(>Chi)'][length(grpsize_icc_ngrps50_lrt[,'Pr(>Chi)'])]
    grpsize_icc_ngrps100_p = grpsize_icc_ngrps100_lrt[,'Pr(>Chi)'][length(grpsize_icc_ngrps100_lrt[,'Pr(>Chi)'])]
    
    #append to output list
    all_int_unpack = list("ngrps_grpsize_icc1" = ngrps_grpsize_icc1_p,
                          "ngrps_grpsize_icc2" = ngrps_grpsize_icc2_p,
                          "ngrps_grpsize_icc3" = ngrps_grpsize_icc3_p,
                          "ngrps_icc_grpsize5" = ngrps_icc_grpsize5_p,
                          "ngrps_icc_grpsize30" = ngrps_icc_grpsize30_p, 
                          "ngrps_icc_grpsize50" = ngrps_icc_grpsize50_p,
                          "grpsize_icc_ngrps30" = grpsize_icc_ngrps30_p,
                          "grpsize_icc_ngrps50" = grpsize_icc_ngrps50_p,
                          "grpsize_icc_ngrps100" = grpsize_icc_ngrps100_p)
    
    outList[["all_int_unpacked"]] = all_int_unpack
    
  }
  
  #Now we're going to test all possible two way interactions (this is assuming the highest level guy isn't significant)
  #get their p-values and do the same procedure as above
  p_val_all_int_grpsize_ngrps = all_int_lrt['ngrps:grpsize', 'Pr(>Chi)']
  p_val_all_int_icc_grpsize = all_int_lrt['icc:grpsize', 'Pr(>Chi)']
  p_val_all_int_icc_ngrps = all_int_lrt['icc:ngrps', 'Pr(>Chi)']

  
  #start with grpsize:ngrps
  if (p_val_all_int_grpsize_ngrps >= .01) {
    
    int_grpsize_ngrps_test_msg = "not significant"
    outList[["grpsize_ngrps_result"]] = int_grpsize_ngrps_test_msg
    
  } else {
    
    int_grpsize_ngrps_test_msg = sprintf('significant, p = %s', p_val_all_int_grpsize_ngrps)
    outList[["grpsize_ngrps_result"]] = int_grpsize_ngrps_test_msg
    
    #test ngrps at various ICC levels
    ngrps_grpsize5 = glm(non_coverage ~ ngrps, data = sim_datasets[sim_datasets$grpsize == "grpsize5",], family = "binomial")
    ngrps_grpsize30 = glm(non_coverage ~ ngrps, data = sim_datasets[sim_datasets$grpsize == "grpsize30",], family = "binomial")
    ngrps_grpsize50 = glm(non_coverage ~ ngrps, data = sim_datasets[sim_datasets$grpsize == "grpsize50",], family = "binomial")
    
    #test ICC at various ngrps levels
    grpsize_ngrps30 = glm(non_coverage ~ grpsize, data = sim_datasets[sim_datasets$ngrps == "ngrps30",], family = "binomial")
    grpsize_ngrps50 = glm(non_coverage ~ grpsize, data = sim_datasets[sim_datasets$ngrps == "ngrps50",], family = "binomial")
    grpsize_ngrps100 = glm(non_coverage ~ grpsize, data = sim_datasets[sim_datasets$ngrps == "ngrps100",], family = "binomial")
    
    #compute likelihood ratio tests
    ngrps_grpsize5_lrt = anova(ngrps_grpsize5, test = "LRT")
    ngrps_grpsize30_lrt = anova(ngrps_grpsize30, test = "LRT")
    ngrps_grpsize50_lrt = anova(ngrps_grpsize50, test = "LRT")
    
    grpsize_ngrps30_lrt = anova(grpsize_ngrps30, test = "LRT")
    grpsize_ngrps50_lrt = anova(grpsize_ngrps50, test = "LRT")
    grpsize_ngrps100_lrt  = anova(grpsize_ngrps100, test = "LRT")
    
    #get pvalues
    ngrps_grpsize5_p = ngrps_grpsize5_lrt[,'Pr(>Chi)'][2]
    ngrps_grpsize30_p = ngrps_grpsize30_lrt[,'Pr(>Chi)'][2]
    ngrps_grpsize50_p = ngrps_grpsize50_lrt[,'Pr(>Chi)'][2]
    
    grpsize_ngrps30_p = grpsize_ngrps30_lrt[,'Pr(>Chi)'][2]
    grpsize_ngrps50_p = grpsize_ngrps50_lrt[,'Pr(>Chi)'][2]
    grpsize_ngrps100_p = grpsize_ngrps100_lrt[,'Pr(>Chi)'][2]
    
    #append to output list
    int_grpsize_ngrps_unpack = list("ngrps_grpsize5" = ngrps_grpsize5_p,
                                    "ngrps_grpsize30" = ngrps_grpsize30_p,
                                    "ngrps_grpsize50" = ngrps_grpsize50_p,
                                    "grpsize_ngrps30" = grpsize_ngrps30_p,
                                    "grpsize_ngrps50" = grpsize_ngrps50_p, 
                                    "grpsize_ngrps100" = grpsize_ngrps100_p)
    
    outList[["grpsize_ngrps_unpack"]] = int_grpsize_ngrps_unpack
    
  }
  
  #move on to icc:grpsize
  if (p_val_all_int_icc_grpsize >= .01) {
    
    int_icc_grpsize_test_msg = "not significant"
    outList[["icc_grpsize_result"]] = int_icc_grpsize_test_msg
    
  } else {
    
    int_icc_grpsize_test_msg = sprintf('significant, p = %s', p_val_all_int_icc_grpsize)
    outList[["icc_grpsize_result"]] = int_icc_grpsize_test_msg
    
    #test grpsize at various ICC levels
    grpsize_icc1 = glm(non_coverage ~ grpsize, data = sim_datasets[sim_datasets$icc == 0.1,], family = "binomial")
    grpsize_icc2 = glm(non_coverage ~ grpsize, data = sim_datasets[sim_datasets$icc == 0.2,], family = "binomial")
    grpsize_icc3 = glm(non_coverage ~ grpsize, data = sim_datasets[sim_datasets$icc == 0.3,], family = "binomial")
    
    #test ICC at various grpsize levels
    icc_grpsize5 = glm(non_coverage ~ icc, data = sim_datasets[sim_datasets$grpsize == "grpsize5",], family = "binomial")
    icc_grpsize30 = glm(non_coverage ~ icc, data = sim_datasets[sim_datasets$grpsize == "grpsize30",], family = "binomial")
    icc_grpsize50 = glm(non_coverage ~ icc, data = sim_datasets[sim_datasets$grpsize == "grpsize50",], family = "binomial")
    
    #compute likelihood ratio tests
    grpsize_icc1_lrt = anova(grpsize_icc1, test = "LRT")
    grpsize_icc2_lrt = anova(grpsize_icc2, test = "LRT")
    grpsize_icc3_lrt = anova(grpsize_icc3, test = "LRT")
    
    icc_grpsize5_lrt = anova(icc_grpsize5, test = "LRT")
    icc_grpsize30_lrt = anova(icc_grpsize30, test = "LRT")
    icc_grpsize50_lrt = anova(icc_grpsize50, test = "LRT")
    
    #get pvalues
    grpsize_icc1_p = grpsize_icc1_lrt[,'Pr(>Chi)'][2]
    grpsize_icc2_p = grpsize_icc2_lrt[,'Pr(>Chi)'][2]
    grpsize_icc3_p = grpsize_icc3_lrt[,'Pr(>Chi)'][2]
    
    icc_grpsize5_p = icc_grpsize5_lrt[,'Pr(>Chi)'][2]
    icc_grpsize30_p = icc_grpsize30_lrt[,'Pr(>Chi)'][2]
    icc_grpsize50_p = icc_grpsize50_lrt[,'Pr(>Chi)'][2]
    
    #append to output list
    int_icc_grpsize_unpack = list("icc_grpsize5" = icc_grpsize5_p,
                                  "icc_grpsize30" = icc_grpsize30_p,
                                  "icc_grpsize50" = icc_grpsize50_p,
                                  "grpsize_icc1" = grpsize_icc1_p,
                                  "grpsize_icc2" = grpsize_icc2_p, 
                                  "grpsize_icc3" = grpsize_icc3_p)
    
    outList[["icc_grpsize_unpack"]] = int_icc_grpsize_unpack
    
  }
  
  #end with icc:ngrps
  if (p_val_all_int_icc_ngrps >= .01) {
    
    int_icc_ngrps_test_msg = "not significant"
    outList[["icc_ngrps_result"]] = int_icc_ngrps_test_msg
    
  } else {
    
    int_icc_ngrps_test_msg = sprintf('significant, p = %s', p_val_all_int_icc_ngrps)
    outList[["icc_ngrps_result"]] = int_icc_ngrps_test_msg
    
    #test ngrps at various ICC levels
    ngrps_icc1 = glm(non_coverage ~ ngrps, data = sim_datasets[sim_datasets$icc == 0.1,], family = "binomial")
    ngrps_icc2 = glm(non_coverage ~ ngrps, data = sim_datasets[sim_datasets$icc == 0.2,], family = "binomial")
    ngrps_icc3 = glm(non_coverage ~ ngrps, data = sim_datasets[sim_datasets$icc == 0.3,], family = "binomial")
    
    #test ICC at various ngrps levels
    icc_ngrps30 = glm(non_coverage ~ icc, data = sim_datasets[sim_datasets$ngrps == "ngrps30",], family = "binomial")
    icc_ngrps50 = glm(non_coverage ~ icc, data = sim_datasets[sim_datasets$ngrps == "ngrps50",], family = "binomial")
    icc_ngrps100 = glm(non_coverage ~ icc, data = sim_datasets[sim_datasets$ngrps == "ngrps100",], family = "binomial")
    
    #compute likelihood ratio tests
    ngrps_icc1_lrt = anova(ngrps_icc1, test = "LRT")
    ngrps_icc2_lrt = anova(ngrps_icc2, test = "LRT")
    ngrps_icc3_lrt = anova(ngrps_icc3, test = "LRT")
    
    icc_ngrps30_lrt = anova(icc_ngrps30, test = "LRT")
    icc_ngrps50_lrt = anova(icc_ngrps50, test = "LRT")
    icc_ngrps100_lrt = anova(icc_ngrps100, test = "LRT")
    
    #get pvalues
    ngrps_icc1_p = ngrps_icc1_lrt[,'Pr(>Chi)'][2]
    ngrps_icc2_p = ngrps_icc2_lrt[,'Pr(>Chi)'][2]
    ngrps_icc3_p = ngrps_icc3_lrt[,'Pr(>Chi)'][2]
    
    icc_ngrps30_p = icc_ngrps30_lrt[,'Pr(>Chi)'][2]
    icc_ngrps50_p = icc_ngrps50_lrt[,'Pr(>Chi)'][2]
    icc_ngrps100_p = icc_ngrps100_lrt[,'Pr(>Chi)'][2]
    
    #append to output list
    int_icc_ngrps_unpack = list("icc_ngrps30" = icc_ngrps30_p,
                                "icc_ngrps50" = icc_ngrps50_p,
                                "icc_ngrps100" = icc_ngrps100_p,
                                "ngrps_icc1" = ngrps_icc1_p,
                                "ngrps_icc2" = ngrps_icc2_p, 
                                "ngrps_icc3" = ngrps_icc3_p)
    
    outList[["icc_ngrps_unpack"]] = int_icc_ngrps_unpack
    
  }
  
  #return nested output list with p-values and significance messages
  #goal is to first summarize significance of interactions so user can do follow-up tests on their own

  return(outList)
  
}