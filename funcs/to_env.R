##This function was taken from the mlmpower package (https://github.com/bkeller2/mlmpower)
##See Enders, Keller, & Woller 2023, Psych Methods (doi: https://doi.org/10.1037/met0000614)

to_env <- function(model) {
  
  # Compute with model
  l <- with(model, {
    
    # Obtain level-1 and level-2 vars
    l1 <- which(vapply(predictors, levels, numeric(1L)) == 1)
    l2 <- which(vapply(predictors, levels, numeric(1L)) == 2)
    
    # Selectors for products
    sel_p <- vapply(actions, \(.) .$type == "product", logical(1L))
    # Selectors for random slopes
    sel_r <- vapply(actions, \(.) .$type == "random_slope", logical(1L))
    
    # Make product weight matrix
    # NOTE Assumes first name is always level-1
    pmat <- matrix(0, length(l2), length(l1))
    for (a in actions[sel_p]) {
      j <- which(names(l1) %in% a$name[1])
      i <- which(names(l2) %in% a$name[2])
      pmat[i, j] <- a$weight
    }
    # Make random slope weights
    rvec <- vector("numeric", length(l1))
    for (a in actions[sel_r]) {
      i <- which(names(l1) %in% a$name)
      rvec[i] <- a$weight
    }
    
    # Create list
    within(list(), {
      # Get means
      mean_Y <- outcome$mean
      names(mean_Y) <- outcome$name
      mean_X <- vapply(predictors[l1], \(.) .$mean, numeric(1L))
      mean_Z <- vapply(predictors[l2], \(.) .$mean, numeric(1L))
      
      # Get model implied means of X based on centering
      model_mean_X <- ifelse(
        vapply(predictors[l1], \(.) 'mp_timevar' %in% class(.), logical(1L)),
        mean_X,
        0
      )
      # Get variances
      var_Y <- outcome$sd^2
      var_X <- vapply(predictors[l1], \(.) .$sd^2, numeric(1L))
      var_Z <- vapply(predictors[l2], \(.) .$sd^2, numeric(1L))
      # Get correlations
      corr_X <- corrs$within
      corr_Z <- corrs$between
      corr_X_Z <- corrs$between
      corr_raneffects <- corrs$randeff
      # Get icc (Assumes only one)
      icc_Y <- if (is.null(outcome$icc)) effect_size$icc else outcome$icc
      icc_X <- vapply(predictors[l1], \(.) {
        if (is.null(.$icc)) effect_size$icc else .$icc
      }, numeric(1L))
      # Get R2s
      R2_X_w <- effect_size$within
      R2_increment_b <- effect_size$between
      R2_XXproduct_w <- NULL # Place holders (not used currently)
      R2_ZZproduct_b <- NULL # Place holders (not used currently)
      R2_XZproduct_w <- effect_size$product
      R2_ranslopes_w <- effect_size$random_slope
      
      # Get weights
      weights_X_w <- vapply(predictors[l1], \(.) .$weight, numeric(1L))
      weights_XXproduct_w <- NULL # Place holders (not used currently)
      weights_ZZproduct_b <- NULL # Place holders (not used currently)
      weights_XZproduct_w <- c(pmat)
      weights_ranslopes_w <- rvec
      weights_increment_b <- c(
        vapply(icc_X, \(.) if (. == 0) 0 else NA, numeric(1L)),
        vapply(predictors[l2], \(.) .$weight, numeric(1L))
      )
      
      # Set binary Z for level2
      binary_Z <- vapply(predictors[l2], \(.) {
        if ("mp_binary" %in% class(.)) .$mean else 0.0
      }, numeric(1L))
    })
  })
  
  # Return list as env
  return(list2env(l))
}