# function by Nicolas Staedler
logit_mod <- function(x, eps = 0.01) {
  x[x == 1] <- 1 - eps
  x[x == 0] <- 0 + eps
  log(x / (1 - x))
}
