########################################################################
###################### M&H2005 Replication Script ######################
########### JF Guassi Moreira (jguassimoreira[at]ucla[dot]edu) #########
########################### Winter 2021/22 #############################

#########
## Documentation
#########

#This is the main script for replicating M&H's 2005 findings
#the following libraries are called (python style), not loaded:
#dplyr
#lme4
#lmerTest
#MASS
#these packages are used in the helper functions used in this script

#########
## Import helper functions, set path to data
#########
studyPath = file.path("C:", "Users", "Joao Guassi Moreira", "Documents", "mh2005_replisim") #declare the path to the study
scriptPath = file.path(sprintf("%s", studyPath), "funcs") #path to helper funcs
file.sources = Sys.glob(file.path(sprintf("%s", scriptPath), "*.R")) #Construct paths to helper funcs, save in vector
invisible(sapply(file.sources,FUN=source)) #use sapply to source each file path in the helper func vector
datPath = file.path(sprintf("%s", studyPath), "data") #path where we'll eventually save our data

#########
## Define the simulation conditions
#########

#these values can be found on the right hand column of M&H2005 pg 88
icc_conds = c(0.1, 0.2, 0.3) #the ICC levels
ngroup_conds = c(30, 50, 100) #L2 sample size levels (using diff var name to differentiate from function vars)
groupsize_conds = c(5, 30, 50) #L1 sample size levels (same note re naming; this naming is consistent w M&H paper)
sim_conds = expand.grid(icc_conds, ngroup_conds, groupsize_conds); names(sim_conds) = c('icc', 'ngroups', 'groupsize') #create dataframe with all conditions

#########
## Generate the x and z data
#########
set.seed(1234) #setting a seed for reproducibility

#M&H used a fixed dataset for their sims since multilevel models assume explanatory variables are fixed
zdat = data.frame(subID = c(paste0('sub',1:30)), z = rnorm(30)) #L2/between cluster variables
xdat = data.frame(subID = rep(paste0('sub',1:30), each = 5), x = rnorm(5*30)) #L1/within cluster variables

#########
## Run the simulation
#########

for (i in 1:dim(sim_conds)[1]) { #loop over the simulations
  
  print(sim_conds[i,])#print current condition
  outList = list() #initialize empty list to hold output
  
  for (r in 1:1000) { #1000 replications per condition, as specified by M&H, pg 89
    outList[r] = sim_mlm_legacy(N_l2 = sim_conds[i, 'ngroups'], #simulate data, save output
                 n_l1 = sim_conds[i, 'groupsize'],
                 icc = sim_conds[i, 'icc'],
                 sigma2 = 1,
                 xdat = xdat, 
                 zdat = zdat)
    
    if((r / 100) == round(r / 100)) {print(sprintf('completed iteration %s', r))} #notify us of the sims progress every hundred reps
    
  }
  
  saveRDS(outList, file = sprint('%s/icc%s_ngrps%s_grpsize%s', datPath, sim_conds[i,'icc']*10, sim_conds[i,'ngroups'], sim_conds[i,'groupsize']))
  
}
