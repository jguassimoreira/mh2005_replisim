########################################################################
################# M&H2005 Replication Analysis Script ##################
########### JF Guassi Moreira (jguassimoreira[at]ucla[dot]edu) #########
########################### Winter 2021/22 #############################

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
