########################################################################
###################### M&H2005 Replication Script ######################
########### JF Guassi Moreira (jguassimoreira[at]ucla[dot]edu) #########
######################### Winter/Spring 2024 ###########################

library(mlmpower) #would prefer to call this like a pythong library but it's much easier to load up front here

#########
## Import helper functions, set path to data
#########
studyPath = file.path("/Users", "jguassimoreira", "Documents", "mh2005_replisim") #declare the path to the study
scriptPath = file.path(sprintf("%s", studyPath), "funcs") #path to helper funcs
file.sources = Sys.glob(file.path(sprintf("%s", scriptPath), "*.R")) #Construct paths to helper funcs, save in vector
invisible(sapply(file.sources,FUN=source)) #use sapply to source each file path in the helper func vector
datPath = file.path(sprintf("%s", studyPath), "data", "extension") #path where we'll eventually save our data
vcov.default = function(object,...) object$vcov #temp

#########
## Define the simulation conditions
#########

#these values are a combination of the original M&H2005 sim (right hand column of M&H2005 pg 88)
#and our proposed extension
icc_conds = c(0.1, 0.2, 0.3, 0.4) #the ICC levels - 0.4 is an extension
ngroup_conds = c(30, 50, 100) #L2 sample size levels (using diff var name to differentiate from function vars)
groupsize_conds = c(5, 30, 50) #L1 sample size levels (same note re naming; this naming is consistent w M&H paper)
imbalance_conds = c(1, 1.5, 2, 2.5) #imbalance between largest and smallest group sizes (extension)
rfx_corr_conds = c(-.6, -.3, 0, .3, .6) #correlation between random slope and intercepts (extension)
mod_conds = c(1, 2, 3) #new model conditions; 1 is original, others add m additional L1 rfx slopes and cross-lev ints
sim_conds = expand.grid(icc_conds, #create dataframe with all conditions
                        ngroup_conds,
                        groupsize_conds,
                        imbalance_conds,
                        rfx_corr_conds,
                        mod_conds) 
names(sim_conds) = c('icc', 'ngroups', 'groupsize', 'imbalance', 'rfx_corr', 'mod') #add column names
sim_conds[sim_conds$mod > 1, 'rfx_corr'] = 0 #set rfx_corr to zero for cells that have more than 1 L1 slope
#the step above ends up creating a bunch of duplicate conditions...
sim_conds = dplyr::distinct(sim_conds) #...and this line gets rid of them


#########
## Create list of mlmpower objects
#########

#coefficient values will be set, and data will be generated, differently from M&H 2005
#first, coefficients will be generated as a function of the percent variance explained by L1, L2 and cross-level
#interactions (Rights & Sterba, 2019). Data will then be generated as specified by Enders et al., 2023, using some of
#the mlmpower package's functionalities. 

#this part creates mlm power object for each condition in the simulation, puts them in a list, and saves them
#I'm grouping model objects into three lists, each based on the propotion variance explained by the predictors:
#high (Level-1 - 0.095, Level-2 - 0.095, cross-level - 0.095) 
#mid (Level-1 - 0.065, Level-2 - 0.065, cross-level - 0.025)
#low (Level-1 - 0.035, Level-2 - 0.035, cross-level - 0.01)
#the high grouping corresponds to the original M&H sim

#NOTES: for model configurations that have more than 1 within var and cross-level (e.g., mod_conds == 2 or 3), (i) raneff 
#cor is set to zero for simplicity and efficiency and (ii) the contribution of each var or int to variance explained is equal

#high variance explained
model_objects_high_var_exp = lapply(1:dim(sim_conds)[1],
                                    FUN = function(x) {make_mlmpwr_objects(mod = sim_conds[x,'mod'], #model configuration
                                                                           sim_icc = sim_conds[x,'icc'], #Y var ICC
                                                                           w_r2 = 0.095, #var explained by within var(s)
                                                                           b_r2 = 0.095, #var explained by btwn var(s)
                                                                           cl_r2 = 0.095, #var explained by cross-lev int(s)
                                                                           raneff_cor = sim_conds[x,'rfx_corr'])}) #correlation btwn random slope & int

#mid variance explained
model_objects_mid_var_exp = lapply(1:dim(sim_conds)[1],
                                    FUN = function(x) {make_mlmpwr_objects(mod = sim_conds[x,'mod'], #model configuration
                                                                           sim_icc = sim_conds[x,'icc'], #Y var ICC
                                                                           w_r2 = 0.065, #var explained by within var(s)
                                                                           b_r2 = 0.065, #var explained by btwn var(s)
                                                                           cl_r2 = 0.025, #var explained by cross-lev int(s)
                                                                           raneff_cor = sim_conds[x,'rfx_corr'])}) #correlation btwn random slope & int

#low variance explained
model_objects_low_var_exp = lapply(1:dim(sim_conds)[1],
                                    FUN = function(x) {make_mlmpwr_objects(mod = sim_conds[x,'mod'], #model configuration
                                                                           sim_icc = sim_conds[x,'icc'], #Y var ICC
                                                                           w_r2 = 0.035, #var explained by within var(s)
                                                                           b_r2 = 0.035, #var explained by btwn var(s)
                                                                           cl_r2 = 0.01, #var explained by cross-lev int(s)
                                                                           raneff_cor = sim_conds[x,'rfx_corr'])}) #correlation btwn random slope & int

#This code runs quickly, so it's not imperative to save these lists to cut down on time, but I think it helps
#for data sharing
saveRDS(model_objects_high_var_exp, file = file.path(sprintf("%s", datPath), "model_objects/high_var_exp.rds"))
saveRDS(model_objects_mid_var_exp, file = file.path(sprintf("%s", datPath), "model_objects/mid_var_exp.rds"))
saveRDS(model_objects_low_var_exp, file = file.path(sprintf("%s", datPath), "model_objects/low_var_exp.rds"))

#########
## Run simulation
#########

#Note that unlike the replication of the original ('main_replication.R'), data generation will occur inside the sim_mlm_extension.R
#function. Moreover, data are simulated once for every rep, unlike the original sim

#start with the high variance explained models
start.time <- Sys.time()
for (i in 1:dim(sim_conds)[1]) {
  
  print(sim_conds[i,])#print current condition
  outList = list() #initialize empty list to hold output
  
  for (r in 1:10) {
    
    outList[[r]] = sim_mlm_extension(N_l2 = sim_conds[i,'ngroups'],
                      n_l1 = sim_conds[i,'groupsize'], 
                      mod_obj = model_objects_high_var_exp[[r]],
                      imbalance = sim_conds[i,'imbalance'])
    names(outList[[r]]) = c('fx_coefs', 'varComp', 'confInt')
    
  }
  
}
end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken


### TEMP STUFF BELOW THAT I NEED TO REMOVE LATER

n_within = 50
n_between = 100
ndata = 1 #don't need -- check that it works with more complex model configurations

lvls <- vapply(mh$predictors, levels, numeric(1L)) #don't need

#don't need
timevar_l1 <- unlist(lapply(
  mh$predictors[lvls == 1],  # Subset level-1
  \(.)  'mp_timevar' %in% class(.)  # Select timevar
))

#don't need
if (TRUE %in% timevar_l1) {
  len <- length(mh$predictors[lvls == 1][[which(timevar_l1)]]$values)
  if (missing(n_within)) {
    n_within <- len
  } else if (len != n_within) {
    cli::cli_alert_info(
      "Setting n_within = {len} to match time variable's length"
    )
    n_within <- len
  }
}

#don't need
if (ndata > 1) {
  return(
    lapply(seq_len(ndata), \(.) {
      model |> generate(n_within, n_between, ndata = 1, mechanism)
    })
  )
}

#don't need
binary_l2 <- sapply(
  # Subset level-2
  mh$predictors[lvls == 2],  # Subset level-2
  \(.)  'mp_binary' %in% class(.)   # Select binary
)

#need
p = make_parameters(mh)

#need
N <- n_within * n_between # Total sample size
l1 <- length(p$mean_X)
l2 <- length(p$mean_Z)

# Generate ID variable
`_id` <- seq_len(n_between) %x% matrix(1, n_within)

# Generate X variable
X_w <- rmvnorm_nomean(N, p$phi_w) #x - within
X_b <- rmvnorm(N, c(p$mean_X, rep(0, l2)), p$phi_b)[`_id`, , drop = F] #x - between


# Create all interactions
inter <- do.call('cbind', apply(
  X_w, 2,
  \(.). * X_b[, seq_len(l2) + l1, drop = F],
  simplify = F
))

# Create predictors matrix
X <- cbind(1, X_w, inter, X_b)

#TO DO -- HERE IS WHERE WE WANT TO ADD IN THE IMBALANCE; I.E., SIMULATE DATA IN A WAY THAT REFLECTS EACH CASE HAS A VALUE AND THEN 'DELETE' SOME RANDOMLY TO INTRODUCE IMBALANCE

# Generate level-1 residuals
e_ij <- rnorm(N, 0, sqrt(p$var_e))

# Generate level-2 residuals
u_j <- rmvnorm_nomean(n_between, p$tau)[`_id`, , drop = F]

# Generate Outcome
Y <- X %*% p$gammas + rowSums(X[,seq_len(l1 + 1), drop = F] * u_j) + e_ij

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

names(d) <- c(
  '_id',
  names(p$mean_Y),
  names(which(vapply(mh$predictors, levels, numeric(1L)) == 1)),
  names(which(vapply(mh$predictors, levels, numeric(1L)) == 2))
)
