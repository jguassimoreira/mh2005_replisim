########################################################################
################# M&H2005 Replication Analysis Script ##################
########### JF Guassi Moreira (jguassimoreira[at]ucla[dot]edu) #########
############################ Spring 2024 ###############################

#########
## Libraries used (others are called with '::')
#########

library(ggplot2)
library(aod)

#########
## Import helper functions, set path to data
#########
studyPath = file.path("/Users", "jguassimoreira", "Documents", "mh2005_replisim") #declare the path to the study
scriptPath = file.path(sprintf("%s", studyPath), "funcs") #path to helper funcs
file.sources = Sys.glob(file.path(sprintf("%s", scriptPath), "*.R")) #Construct paths to helper funcs, save in vector
invisible(sapply(file.sources,FUN=source)) #use sapply to source each file path in the helper func vector
datPath = file.path(sprintf("%s", studyPath), "data") #path where we'll eventually save our data

#########
## Analyze the sim
#########

#read in all the files
conds_to_load = Sys.glob(file.path(sprintf("%s/%s", datPath, "replication"), "*.RDS")) #get file paths for all M&H simulation conditions
conds_list = lapply(conds_to_load, FUN = readRDS) #read in the conditions

params_to_analyze = c("int_var", "x_var", "sigma2", "intercept", "b_x", "b_z", "b_xz") #analyze the parameters

#list of analyses, go through each element to find results
analysisList = lapply(1:length(params_to_analyze), function(x) analyze_sim(conds_list, params_to_analyze[x], conds_to_load))

#now do the interaction
interactionList = lapply(1:length(params_to_analyze), function(x) replication_interaction_analysis(conds_list, params_to_analyze[x], conds_to_load))

## results
#INTERCEPT VARIANCE
#No sig 3 way int

#Icc x grpsize is sig (p = .001047)
#Effect of ICC is sig when grpsize = 5 (p = 1.23 e-08)
#Effect of grpsize is sig when ICC = 0.1 (p = 1.52 e-09)
#Effect of grpsize is sig when ICC = 0.2 (p = .00159)

#X VARIANCE
#3 way int is sig, p = 0.00267
#Ngrps x grpsize is sig @ ICC1 (p = .008)
#Ngrps x grpsize is sig @ ICC2 (p = .002)
#ngrps x icc is sig @ grpsize = 5 (p = .001)
#grpsize x icc is sig at ngrps = 100 (p = 2.703e-05)

#sigma2
#nothing sig

#intercept
#nothing sig

#fixed effect of x
#nothing sig

#fixed effect of z
#nothing sig

#cross level interaction
#nothing sig

### Follow-up on intercept variance results
sim_data_int_var = interactionList[[1]][['data']] 

#We know effect of ICC is sig when grpsize = 5, so let's see what that looks like
aggregate(non_coverage ~ icc, FUN = "mean", data = sim_data_int_var[sim_data_int_var$grpsize == 'grpsize5',])

#Now do effect of grpsize when ICC = 0.1
aggregate(non_coverage ~ grpsize, FUN = "mean", data = sim_data_int_var[sim_data_int_var$icc == '0.1',])

#Now do effect of grpsize when ICC = 0.2
aggregate(non_coverage ~ grpsize, FUN = "mean", data = sim_data_int_var[sim_data_int_var$icc == '0.2',])

### Follow-up on random slope variance (x variance)
sim_data_x_var = interactionList[[2]][['data']] 

#Ngrps x grpsize is sig @ ICC1 (p = .008)
#test ngrps at grpsize 5 - significant
anova(glm(non_coverage ~ ngrps, data = sim_data_x_var[sim_data_x_var$icc == "0.1" & sim_data_x_var$grpsize == "grpsize5",], family = "binomial"), test = "LRT")
aggregate(non_coverage ~ ngrps, FUN = "mean", data = sim_data_x_var[sim_data_x_var$icc == "0.1" & sim_data_x_var$grpsize == "grpsize5",])
#test ngrps at grpsize 30 - not significant
anova(glm(non_coverage ~ ngrps, data = sim_data_x_var[sim_data_x_var$icc == "0.1" & sim_data_x_var$grpsize == "grpsize30",], family = "binomial"), test = "LRT")
#test ngrps at grpsize 50 - not significant
anova(glm(non_coverage ~ ngrps, data = sim_data_x_var[sim_data_x_var$icc == "0.1" & sim_data_x_var$grpsize == "grpsize50",], family = "binomial"), test = "LRT")

#test grpsize at ngrps 30 - significant
anova(glm(non_coverage ~ grpsize, data = sim_data_x_var[sim_data_x_var$icc == "0.1" & sim_data_x_var$ngrps == "ngrps30",], family = "binomial"), test = "LRT")
aggregate(non_coverage ~ grpsize, FUN = "mean", data = sim_data_x_var[sim_data_x_var$icc == "0.1" & sim_data_x_var$ngrps == "ngrps30",])
#test grpsize at ngrps 50 - significant
anova(glm(non_coverage ~ grpsize, data = sim_data_x_var[sim_data_x_var$icc == "0.1" & sim_data_x_var$ngrps == "ngrps50",], family = "binomial"), test = "LRT")
aggregate(non_coverage ~ grpsize, FUN = "mean", data = sim_data_x_var[sim_data_x_var$icc == "0.1" & sim_data_x_var$ngrps == "ngrps50",])
#test grpsize at ngrps 100 - significant
anova(glm(non_coverage ~ grpsize, data = sim_data_x_var[sim_data_x_var$icc == "0.1" & sim_data_x_var$ngrps == "ngrps100",], family = "binomial"), test = "LRT")
aggregate(non_coverage ~ grpsize, FUN = "mean", data = sim_data_x_var[sim_data_x_var$icc == "0.1" & sim_data_x_var$ngrps == "ngrps100",])

#ngrps x grpsize is sig @ ICC2 (p = .002)
#test grpsize at ngrps 30 - significant
anova(glm(non_coverage ~ grpsize, data = sim_data_x_var[sim_data_x_var$icc == "0.2" & sim_data_x_var$ngrps == "ngrps30",], family = "binomial"), test = "LRT")
aggregate(non_coverage ~ grpsize, FUN = "mean", data = sim_data_x_var[sim_data_x_var$icc == "0.2" & sim_data_x_var$ngrps == "ngrps30",])
#test grpsize at ngrps 50 
anova(glm(non_coverage ~ grpsize, data = sim_data_x_var[sim_data_x_var$icc == "0.2" & sim_data_x_var$ngrps == "ngrps50",], family = "binomial"), test = "LRT")
#test grpsize at ngrps 100 
anova(glm(non_coverage ~ grpsize, data = sim_data_x_var[sim_data_x_var$icc == "0.2" & sim_data_x_var$ngrps == "ngrps100",], family = "binomial"), test = "LRT")

#test ngrps at grpsize 5
anova(glm(non_coverage ~ ngrps, data = sim_data_x_var[sim_data_x_var$icc == "0.2" & sim_data_x_var$grpsize == "grpsize5",], family = "binomial"), test = "LRT")
aggregate(non_coverage ~ ngrps, FUN = "mean", data = sim_data_x_var[sim_data_x_var$icc == "0.2" & sim_data_x_var$grpsize == "grpsize5",])
#test ngrps at grpsize 30 - not significant
anova(glm(non_coverage ~ ngrps, data = sim_data_x_var[sim_data_x_var$icc == "0.2" & sim_data_x_var$grpsize == "grpsize30",], family = "binomial"), test = "LRT")
#test ngrps at grpsize 50 - not significant
anova(glm(non_coverage ~ ngrps, data = sim_data_x_var[sim_data_x_var$icc == "0.2" & sim_data_x_var$grpsize == "grpsize50",], family = "binomial"), test = "LRT")

#ngrps x icc is sig @ grpsize = 5 (p = .001)
#test ngrps at icc 0.1 
anova(glm(non_coverage ~ ngrps, data = sim_data_x_var[sim_data_x_var$icc == "0.1" & sim_data_x_var$grpsize == "grpsize5",], family = "binomial"), test = "LRT")
aggregate(non_coverage ~ ngrps, FUN = "mean", data = sim_data_x_var[sim_data_x_var$icc == "0.1" & sim_data_x_var$grpsize == "grpsize5",])
#test ngrps at icc 0.2 
anova(glm(non_coverage ~ ngrps, data = sim_data_x_var[sim_data_x_var$icc == "0.2" & sim_data_x_var$grpsize == "grpsize5",], family = "binomial"), test = "LRT")
aggregate(non_coverage ~ ngrps, FUN = "mean", data = sim_data_x_var[sim_data_x_var$icc == "0.2" & sim_data_x_var$grpsize == "grpsize5",])
#test ngrps at icc 0.3 - not significant
anova(glm(non_coverage ~ ngrps, data = sim_data_x_var[sim_data_x_var$icc == "0.3" & sim_data_x_var$grpsize == "grpsize5",], family = "binomial"), test = "LRT")

#test icc at ngrps 30
anova(glm(non_coverage ~ icc, data = sim_data_x_var[sim_data_x_var$ngrps == "ngrps30" & sim_data_x_var$grpsize == "grpsize5",], family = "binomial"), test = "LRT")
aggregate(non_coverage ~ icc, FUN = "mean", data = sim_data_x_var[sim_data_x_var$ngrps == "ngrps30" & sim_data_x_var$grpsize == "grpsize5",])
#test icc at ngrps 50
anova(glm(non_coverage ~ icc, data = sim_data_x_var[sim_data_x_var$ngrps == "ngrps50" & sim_data_x_var$grpsize == "grpsize5",], family = "binomial"), test = "LRT")
aggregate(non_coverage ~ icc, FUN = "mean", data = sim_data_x_var[sim_data_x_var$ngrps == "ngrps50" & sim_data_x_var$grpsize == "grpsize5",])
#test icc at ngrps 100 
anova(glm(non_coverage ~ icc, data = sim_data_x_var[sim_data_x_var$ngrps == "ngrps100" & sim_data_x_var$grpsize == "grpsize5",], family = "binomial"), test = "LRT")
aggregate(non_coverage ~ icc, FUN = "mean", data = sim_data_x_var[sim_data_x_var$ngrps == "ngrps100" & sim_data_x_var$grpsize == "grpsize5",])

#grpsize x icc is sig at ngrps = 100 (p = 2.703e-05)
#test grpsize at icc 0.1
anova(glm(non_coverage ~ grpsize, data = sim_data_x_var[sim_data_x_var$icc == "0.1" & sim_data_x_var$ngrps == "ngrps100",], family = "binomial"), test = "LRT")
aggregate(non_coverage ~ grpsize, FUN = "mean", data = sim_data_x_var[sim_data_x_var$icc == "0.1" & sim_data_x_var$ngrps == "ngrps100",])
#test grpsize at icc 0.2 - not significant
anova(glm(non_coverage ~ grpsize, data = sim_data_x_var[sim_data_x_var$icc == "0.2" & sim_data_x_var$ngrps == "ngrps100",], family = "binomial"), test = "LRT")
#test grpsize at icc 0.3 - not significant
anova(glm(non_coverage ~ grpsize, data = sim_data_x_var[sim_data_x_var$icc == "0.3" & sim_data_x_var$ngrps == "ngrps100",], family = "binomial"), test = "LRT")

#test icc at grpsize 5
anova(glm(non_coverage ~ icc, data = sim_data_x_var[sim_data_x_var$grpsize == "grpsize5" & sim_data_x_var$ngrps == "ngrps100",], family = "binomial"), test = "LRT")
aggregate(non_coverage ~ icc, FUN = "mean", data = sim_data_x_var[sim_data_x_var$grpsize == "grpsize5" & sim_data_x_var$ngrps == "ngrps100",])
#test icc at grpsize 30 - not significant
anova(glm(non_coverage ~ icc, data = sim_data_x_var[sim_data_x_var$grpsize == "grpsize30" & sim_data_x_var$ngrps == "ngrps100",], family = "binomial"), test = "LRT")
#test icc at grpsize 50 - not significant
anova(glm(non_coverage ~ icc, data = sim_data_x_var[sim_data_x_var$grpsize == "grpsize50" & sim_data_x_var$ngrps == "ngrps100",], family = "binomial"), test = "LRT")

