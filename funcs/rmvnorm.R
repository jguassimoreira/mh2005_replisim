##This function was taken from the mlmpower package (https://github.com/bkeller2/mlmpower)
##See Enders, Keller, & Woller 2023, Psych Methods (doi: https://doi.org/10.1037/met0000614)

rmvnorm <- function(n, mu, Sigma) {
  p <- NCOL(Sigma)
  if (p == 0) return(matrix(0, n, 0))
  (matrix(mu, nrow = n, ncol = p, byrow = T)
    + matrix(rnorm(p * n), ncol = p)
    %*% with(
      eigen(Sigma, symmetric = TRUE), {
        t(vectors %*% (t(vectors) * sqrt(pmax(values, 0))))
      }
    )
  )
}