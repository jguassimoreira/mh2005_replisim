##This function was taken from the mlmpower package (https://github.com/bkeller2/mlmpower)
##See Enders, Keller, & Woller 2023, Psych Methods (doi: https://doi.org/10.1037/met0000614)

diagonal <- function(x) {
  if (is.matrix(x)) diag(x)
  else diag(x, NROW(x), NROW(x))
}
