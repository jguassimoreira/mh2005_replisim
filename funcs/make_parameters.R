##This function was taken from the mlmpower package (https://github.com/bkeller2/mlmpower)
##See Enders, Keller, & Woller 2023, Psych Methods (doi: https://doi.org/10.1037/met0000614)

make_parameters <- function(model) {
  # Convert model to old specification
  model |> to_env() -> env
  
  # Create list using env
  l <- with(env, {
    # variable counts
    num_X <- length(mean_X)
    num_Z <- length(var_Z)
    
    # set binary variable variance
    var_Z[binary_Z > 0 & binary_Z < 1] <- binary_Z[binary_Z > 0 & binary_Z < 1] * (1 - binary_Z[binary_Z > 0 & binary_Z < 1])
    
    # compute within and between-cluster variances
    icc_X[weights_increment_b[seq_len(num_X)] == 0] <- 0
    var_X_b <- var_X * icc_X
    var_X_w <- var_X - var_X_b
    var_Y_b <- var_Y * icc_Y
    var_Y_w <- var_Y - var_Y_b
    
    # construct the within-cluster correlation and covariance matrix of the level-1 Xs
    cor_XX_w <- matrix(1, nrow = num_X, ncol = num_X)
    cor_XX_w[lower.tri(cor_XX_w)] <- cor_XX_w[upper.tri(cor_XX_w)] <- corr_X(num_X * (num_X - 1) / 2)
    phi_XX_w <- diagonal(sqrt(var_X_w)) %*% cor_XX_w %*% diagonal(sqrt(var_X_w))
    
    # construct the between-cluster covariance matrix of X's level-2 group means
    phi_XX_b <- diagonal(sqrt(var_X_b)) %*% cor_XX_w %*% diagonal(sqrt(var_X_b))
    
    # construct the between-cluster correlation and covariance matrix of the level-2 Zs
    cor_ZZ_b <- matrix(1, nrow = num_Z, ncol = num_Z)
    cor_ZZ_b[lower.tri(cor_ZZ_b)] <- cor_ZZ_b[upper.tri(cor_ZZ_b)] <- corr_Z(num_Z * (num_Z - 1) / 2)
    phi_ZZ_b <- diagonal(sqrt(var_Z)) %*% cor_ZZ_b %*% diagonal(sqrt(var_Z))
    
    # for binary level-2 variables, convert inputted correlations to point-biserial correlations then to covariances
    if (is.null(binary_Z)) {
      binary_Z_ind <- rep(FALSE, length(var_Z))
    } else {
      binary_Z_ind <- binary_Z > 0 & binary_Z < 1
    }
    for (r in seq_len(num_Z)) {
      for (i in seq_len(num_Z)) {
        if (binary_Z_ind[r] == T & i != r) {
          cor_pointbi <- cor_ZZ_b[r, i] / sqrt(binary_Z[r] * (1 - binary_Z[r])) * dnorm(qnorm(binary_Z[r]))
          cor_ZZ_b[r, i] <- cor_ZZ_b[i, r] <- cor_pointbi
          phi_ZZ_b[r, i] <- phi_ZZ_b[i, r] <- cor_ZZ_b[r, i] * sqrt(phi_ZZ_b[r, r] * phi_ZZ_b[i, i])
        }
      }
    }
    
    # construct the between-cluster covariance matrix of the X cluster means and level-2 Zs
    cor_XZ_b <- matrix(corr_X_Z(num_X * num_Z), nrow = num_Z, ncol = num_X)
    phi_XZ_b <- diagonal(sqrt(var_Z)) %*% cor_XZ_b %*% diagonal(sqrt(var_X_b))
    
    # construct the between-cluster covariance matrix
    phi_b <- rbind(cbind(phi_XX_b, t(phi_XZ_b)), cbind(phi_XZ_b, t(phi_ZZ_b)))
    
    # construct the within-cluster covariance matrix
    phi_XZwithXZ_w <- phi_XX_w %x% phi_ZZ_b
    phi_XwithXZ_w <- matrix(0, nrow = NROW(phi_XZwithXZ_w), ncol = NCOL(phi_XX_w))
    phi_w <- rbind(cbind(phi_XX_w, t(phi_XwithXZ_w)), cbind(phi_XwithXZ_w, t(phi_XZwithXZ_w)))
    
    # solve for the within-cluster regression coefficients
    weights_scaled <- 1 / sqrt(diagonal(phi_XX_w)) * weights_X_w
    gamma_X_w <- weights_scaled * c(sqrt((var_Y * R2_X_w) / t(weights_scaled) %*% phi_XX_w %*% weights_scaled))
    
    # solve for the cross-level product coefficients
    weights_scaled <- 1 / sqrt(diagonal(phi_XZwithXZ_w)) * weights_XZproduct_w
    if (sum(weights_XZproduct_w) == 0) {
      gamma_XZ_w <- weights_XZproduct_w
    } else {
      gamma_XZ_w <- weights_scaled * c(sqrt((var_Y * R2_XZproduct_w) / t(weights_scaled) %*% phi_XZwithXZ_w %*% weights_scaled))
    }
    gamma_w <- c(gamma_X_w, gamma_XZ_w)
    
    # compute the within-cluster residual variance
    var_e_w <- (1 - icc_Y - R2_X_w - R2_XZproduct_w - R2_ranslopes_w) * var_Y
    
    # solve for the between-cluster coefficients
    select_weighted <- !is.na(weights_increment_b)
    select_nonweighted <- is.na(weights_increment_b)
    if (sum(select_weighted) != NROW(phi_b)) {
      phi_nonweighted_b <- phi_b[select_nonweighted, select_nonweighted, drop = F]
      phi_weighted_b <- phi_b[select_weighted, select_weighted, drop = F]
      phi_covs_b <- phi_b[select_nonweighted, select_weighted, drop = F]
      resvar_Z_b <- phi_weighted_b - t(solve(phi_nonweighted_b) %*% phi_covs_b) %*% phi_covs_b
      
      # Predefine and only compute for non zero phi_b
      weights_scaled <- numeric(sum(select_weighted))
      sel_weight <- diagonal(resvar_Z_b) != 0 # Select out only non-zero diagonals
      weights_scaled[sel_weight] <- 1 / sqrt(diagonal(resvar_Z_b)[sel_weight]) * weights_increment_b[select_weighted][sel_weight]
      gamma_weighted_b <- weights_scaled * c(sqrt((var_Y * R2_increment_b) / t(weights_scaled) %*% resvar_Z_b %*% weights_scaled))
      gamma_b <- c(gamma_w[seq_len(num_X)], rep(0, num_Z))
      gamma_b[select_weighted] <- gamma_weighted_b
      gamma_nonweighted_b <- gamma_b[select_nonweighted]
    } else {
      # Predefine and only compute for non zero phi_b
      weights_scaled <- numeric(NROW(phi_b))
      sel_weight <- diagonal(phi_b) != 0 # Select out only non-zero diagonals
      weights_scaled[sel_weight] <- 1 / sqrt(diagonal(phi_b)[sel_weight]) * weights_increment_b[sel_weight]
      gamma_b <- weights_scaled * c(sqrt((var_Y * R2_increment_b) / t(weights_scaled) %*% phi_b %*% weights_scaled))
    }
    # Replace NaN with 0
    gamma_b[is.na(gamma_b)] <- 0
    
    # compute the random effect correlation matrix
    cor_raneffects <- diag(nrow = num_X + 1, ncol = num_X + 1)
    cor_sel <- which(weights_ranslopes_w != 0)
    cor_raneffects[1, cor_sel + 1] <- cor_raneffects[cor_sel + 1, 1] <- corr_raneffects(length(cor_sel))
    
    # solve for the random slope variances
    cor_ranslopes <- cor_raneffects[-1, -1, drop = F]
    tau_trace <- (var_Y * R2_ranslopes_w) / sum(diagonal(cor_ranslopes %*% phi_XX_w %*% diagonal(weights_ranslopes_w / diagonal(phi_XX_w))))
    var_ranslopes <- if (is.finite(tau_trace)) weights_ranslopes_w * tau_trace / diagonal(phi_XX_w) else rep(0.0, NROW(weights_ranslopes_w))
    
    # vector of correlations between random intercept and slopes
    cor_is <- cor_raneffects[-1, 1, drop = F]
    
    # covariance matrix of just the random slopes
    tau_ranslopes <- diagonal(sqrt(c(var_ranslopes))) %*% cor_raneffects[-1, -1, drop = F] %*% diagonal(sqrt(c(var_ranslopes)))
    
    # compute random intercept variance
    # explained level-2 variation
    b <- t(gamma_b) %*% phi_b %*% gamma_b
    # random intercept variation due to non-zero level-1 means
    a <- t(model_mean_X) %*% (cor_is * sqrt(var_ranslopes))
    s <- t(model_mean_X) %*% tau_ranslopes %*% model_mean_X
    
    # Precompute portion of tau00
    comp1 <- -4 * b * a^2 + a^4 - 4 * a^2 * s + 4 * a^2 * var_Y_b
    # Check if comp1 is positive
    if (is.na(comp1) | comp1 < 0) {
      throw_error(c(
        'The random intercept variance is negative.',
        'x' = 'This is caused because the effect size specified are impossible.',
        'i' = 'The between effect size is most likely too large compared to the ICC.',
        '>' = 'ICC = {icc_Y} and Between R2 = {R2_increment_b}'
      ))
    }
    tau00 <- 0.5 * (-2 * b + a^2 - 2 * s + 2 * var_Y_b) - 0.5 * sqrt(comp1)
    # Check if tau00 is positive
    if (is.na(tau00) | tau00 < 0) {
      throw_error(c(
        'The random intercept variance is negative.',
        'x' = 'This is caused because the effect size specified are impossible.',
        'i' = 'The between effect size is most likely too large compared to the ICC.',
        '>' = 'ICC = {icc_Y} and Between R2 = {R2_increment_b}'
      ))
    }
    
    # compute intercept-slope covariance and construct tau matrix
    tau <- diagonal(sqrt(c(tau00, var_ranslopes))) %*% cor_raneffects %*% diagonal(sqrt(c(tau00, var_ranslopes)))
    cor_raneffects[tau == 0] <- 0
    
    # compute fixed intercept and construct coefficient matrix
    # the mean of the product from Bohrnstedt & Goldberger Equation 3 simplifies because cov(X_w,Z) = 0
    means <- c(model_mean_X, rep(0, num_X + num_Z + num_X*num_Z)) # Set all Z means to 0
    gamma00 <- mean_Y - c(gamma_w, gamma_b) %*% means
    gammas <- c(gamma00, gamma_w, gamma_b)
    
    # variable names
    if (length(names(mean_X)) != 0) {
      vars_Xw <- paste0(names(mean_X), "_w")
      vars_Xb <- paste0(names(mean_X), "_b")
    } else {
      vars_Xw <- vars_Xb <- names(mean_X)
    }
    vars_Z <- names(mean_Z)
    vars_Z[binary_Z_ind == T] <- paste0(vars_Z[binary_Z_ind == T], "_binary")
    vars_XZ <- unlist(lapply(vars_Xw, \(x) {
      vapply(vars_Z, \(w) {
        paste0(x, "*", w)
      }, character(1L))
    }))
    
    # collect parameters and construct names
    params_coeff <- matrix(c(gammas), ncol = 1)
    params_res <- matrix(c(var_e_w), ncol = 1)
    row.names(tau) <- colnames(tau) <- c("Icept", vars_Xw)
    row.names(params_coeff) <- c("Icept", vars_Xw, vars_XZ, vars_Xb, vars_Z)
    row.names(params_res) <- c("Res. Var.")
    colnames(params_coeff) <- colnames(params_res) <- "Value"
    row.names(phi_XX_w) <- colnames(phi_XX_w) <- vars_Xw
    row.names(phi_b) <- colnames(phi_b) <- c(vars_Xb, vars_Z)
    row.names(phi_XZwithXZ_w) <- colnames(phi_XZwithXZ_w) <- vars_XZ
    
    # R-square summary
    check_var_Y <- (
      t(gamma_w) %*% phi_w %*% gamma_w + t(gamma_b) %*%
        phi_b %*% gamma_b + sum(diagonal(tau[-1, -1, drop = F] %*% phi_XX_w))
      + t(model_mean_X) %*% tau_ranslopes %*% model_mean_X + t(model_mean_X) %*%
        tau[-1, 1, drop = F] + tau00 + var_e_w
    )
    # check_var_Y_w <- (
    #     t(gamma_w) %*% phi_w %*% gamma_w
    #     + sum(diag(tau[-1, -1, drop = F] %*% phi_XX_w)) + var_e_w
    # )
    # check_var_Y_b <- (
    #     t(gamma_b) %*% phi_b %*% gamma_b + t(model_mean_X) %*%
    #         tau_ranslopes %*% model_mean_X + t(model_mean_X) %*% tau[-1,1]
    #     + tau00
    # )
    R2check_X_w <- t(gamma_X_w) %*% phi_XX_w %*% gamma_X_w / check_var_Y
    R2check_XZ_w <- t(gamma_XZ_w) %*% phi_XZwithXZ_w %*% gamma_XZ_w / check_var_Y
    R2check_ranslopes_w <- sum(diagonal(tau_ranslopes %*% phi_XX_w)) / check_var_Y
    R2check_var_e <- ((1 - icc_Y) * var_Y - t(gamma_w) %*% phi_w %*% gamma_w - sum(diagonal(tau_ranslopes %*% phi_XX_w))) / check_var_Y
    R2check_XZ_b <- t(gamma_b) %*% phi_b %*% gamma_b / check_var_Y
    if (sum(select_weighted) != NROW(phi_b)) {
      R2check_increment_b <- t(gamma_weighted_b) %*% resvar_Z_b %*% gamma_weighted_b / check_var_Y
    } else {
      R2check_increment_b <- R2check_XZ_b
    }
    R2check_totalminusincrement_b <- R2check_XZ_b - R2check_increment_b
    R2check_ranicept <- tau00 / check_var_Y
    
    # Collect r2 summaries
    r2  <- matrix(
      c(
        R2check_X_w,
        R2check_XZ_w,
        R2check_ranslopes_w,
        R2check_var_e,
        R2check_XZ_b,
        R2check_increment_b,
        R2check_totalminusincrement_b,
        R2check_ranicept
      ),
      dimnames = list(c(
        "Variance Within-Cluster Fixed Effects",
        "Variance Cross-Level Interaction Effects",
        "Variance Random Slopes",
        "Within-Cluster Error Variance",
        "Variance Between-Cluster Fixed Effects",
        "Incremental Variance Level-2 Predictors",
        "Between Variance Level-1 Covariates",
        "Variance Random Intercepts"
      ), 'Proportion')
    )
    
    # Subset phi_XZwithXZ_w based on actual requested products
    phi_XZwithXZ_w <- phi_XZwithXZ_w[gamma_XZ_w != 0, gamma_XZ_w != 0, drop = F]
    
    # Return final output
    list(
      mean_Y = mean_Y,
      gammas = params_coeff,
      tau = tau,
      var_e = params_res,
      mean_X = mean_X,
      mean_Z = mean_Z,
      phi_w = phi_XX_w,
      phi_p = phi_XZwithXZ_w,
      phi_b = phi_b,
      r2 = r2
    )
  })
  
  # Return mp_parameters class
  structure(list2env(l), class = c("mp_parameters", "mp_base"))
}