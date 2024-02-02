########################################################################
###################### M&H2005 Replication Script ######################
########### JF Guassi Moreira (jguassimoreira[at]ucla[dot]edu) #########
######################### Winter/Spring 2024 ###########################

library(mlmpower)

#########
## Import helper functions, set path to data
#########
studyPath = file.path("Users", "jguassimoreira", "Documents", "mh2005_replisim") #declare the path to the study
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
imbalance_conds = c(1, 1.5, 2, 2.5) #Imbalance between largest and smallest group sizes (extension)
sim_conds = expand.grid(icc_conds, ngroup_conds, groupsize_conds); names(sim_conds) = c('icc', 'ngroups', 'groupsize') #create dataframe with all conditions


example1 <- (
  effect_size(
    icc = c(.10,.25),
    within = .065,
    between = .065,
    product = .01,
    random_slope = .03
  )
  + outcome('y', mean = 50, sd = 10)
  + within_predictor('x1', icc = 0, weight = .80)
  + within_predictor('x2', weight = .20)
  + between_predictor('z1', weight = .80)
  + between_predictor('z2', weight = .20)
  + product('x1','z1', weight = 1)
  + random_slope('x1', weight = 1)
)


mh <- (
  effect_size(
    icc = 0.1,
    within = .095,
    between = .095,
    product = .095,
    random_slope = .1
  )
  + outcome('y', mean = 1, sd = .95)
  + within_predictor('x1', icc = 0)
  + between_predictor('z1')
  + product('x1','z1', weight = 1)
  + random_slope('x1', weight = 1)
  + correlations(randeff = 0)
)

n_within = 50
n_between = 100
ndata = 1

lvls <- vapply(mh$predictors, levels, numeric(1L))

timevar_l1 <- unlist(lapply(
  mh$predictors[lvls == 1],  # Subset level-1
  \(.)  'mp_timevar' %in% class(.)  # Select timevar
))

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

if (ndata > 1) {
  return(
    lapply(seq_len(ndata), \(.) {
      model |> generate(n_within, n_between, ndata = 1, mechanism)
    })
  )
}

binary_l2 <- sapply(
  # Subset level-2
  mh$predictors[lvls == 2],  # Subset level-2
  \(.)  'mp_binary' %in% class(.)   # Select binary
)

p = make_parameters(mh)

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
