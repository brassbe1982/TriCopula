R_true <- matrix(c(1, .5, .3, .5, 1, .4, .3, .4, 1), 3)

sim_mixed <- function(n, copula = "gaussian", df = 5, seed = 1) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n))
  U <- tc_rcopula(n, R_true, copula, df)
  d$y1 <- qnorm(U[, 1], 1 + 0.5 * d$x, 1.2)
  d$y2 <- qlnorm(U[, 2], 0.2 + 0.3 * d$x, 0.5)
  d$y3 <- qgamma(U[, 3], shape = 2.5, scale = exp(0.1 + 0.2 * d$x) / 2.5)
  d
}

true_par <- function(copula = "gaussian", df = 5) {
  th <- c(1, 0.5, log(1.2), 0.2, 0.3, log(0.5), 0.1, 0.2, log(2.5),
          TriCopula:::.tc_dep_inverse(R_true))
  if (copula == "student") th <- c(th, log(df - 2.05))
  th
}

fit_mixed <- function(d, copula = "gaussian", ...) {
  tricopula(list(y1 ~ x, y2 ~ x, y3 ~ x), d,
            margins = c("gaussian", "lognormal", "gamma"), copula = copula, ...)
}
