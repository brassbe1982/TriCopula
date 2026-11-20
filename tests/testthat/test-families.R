sim_new <- function(n, copula = "gaussian", df = 5, seed = 1) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n))
  U <- tc_rcopula(n, R_true, copula, df)
  d$y1 <- qweibull(U[, 1], shape = 1.8, scale = exp(0.5 + 0.4 * d$x))
  d$y2 <- TriCopula:::.tc_qinvgauss(U[, 2], exp(0.2 + 0.3 * d$x), 0.6)
  mu <- plogis(-0.3 + 0.5 * d$x)
  d$y3 <- qbeta(U[, 3], mu * 8, (1 - mu) * 8)
  d
}
new_margins <- c("weibull", "invgauss", "beta")
true_new <- function(copula = "gaussian", df = 5) {
  th <- c(0.5, 0.4, log(1.8), 0.2, 0.3, log(0.6), -0.3, 0.5, log(8),
          TriCopula:::.tc_dep_inverse(R_true))
  if (copula == "student") th <- c(th, log(df - 2.05))
  th
}

test_that("new margin densities, CDFs and quantiles are mutually consistent", {
  y <- c(0.2, 0.7, 1.5, 3)
  X <- matrix(1, 4, 1)
  # Weibull against stats
  m <- TriCopula:::.tc_margin_eval(y, X, "weibull", 0.3, log(1.7))
  expect_equal(m$logf, dweibull(y, 1.7, exp(0.3), log = TRUE))
  expect_equal(m$u, pweibull(y, 1.7, exp(0.3)))
  # inverse Gaussian: density integrates to the CDF, quantile inverts it
  m <- TriCopula:::.tc_margin_eval(y, X, "invgauss", 0.3, log(0.6))
  dens <- function(v) exp(TriCopula:::.tc_margin_eval(v, matrix(1, length(v), 1),
                                                      "invgauss", 0.3, log(0.6))$logf)
  for (i in seq_along(y)) {
    expect_equal(integrate(dens, 0, y[i], rel.tol = 1e-10)$value, m$u[i], tolerance = 1e-7)
  }
  expect_equal(integrate(function(v) v * dens(v), 0, Inf)$value, exp(0.3), tolerance = 1e-6)
  expect_equal(integrate(function(v) (v - exp(0.3))^2 * dens(v), 0, Inf)$value,
               exp(0.3)^3 * 0.6^2, tolerance = 1e-5)
  expect_equal(TriCopula:::.tc_qinvgauss(m$u, exp(0.3), 0.6), y, tolerance = 1e-8)
  # Beta against stats
  yb <- c(0.1, 0.4, 0.8, 0.95)
  m <- TriCopula:::.tc_margin_eval(yb, X, "beta", 0.4, log(6))
  mu <- plogis(0.4)
  expect_equal(m$logf, dbeta(yb, mu * 6, (1 - mu) * 6, log = TRUE))
  expect_equal(m$u, pbeta(yb, mu * 6, (1 - mu) * 6))
  # means
  expect_equal(TriCopula:::.tc_margin_mean(0.3, log(1.7), "weibull"),
               exp(0.3) * gamma(1 + 1 / 1.7))
  expect_equal(TriCopula:::.tc_margin_mean(0.4, log(6), "beta"), mu)
})

test_that("margin scores match numerical differentiation for every family", {
  set.seed(2)
  n <- 30
  X <- cbind(1, rnorm(n))
  # two-parameter margins; the three-parameter ones are in test-margins3.R
  for (mg in setdiff(TriCopula:::.tc_margin_names, TriCopula:::.tc_margin_names3)) {
    y <- switch(mg, gaussian = rnorm(n), beta = runif(n, 0.05, 0.95), rexp(n) + 0.1)
    b <- c(0.2, 0.1); e <- 0.3
    f <- function(th, what) {
      m <- TriCopula:::.tc_margin_eval(y, X, mg, th[1:2], th[3])
      if (what == "l") sum(m$logf) else sum(m$u)
    }
    m <- TriCopula:::.tc_margin_eval(y, X, mg, b, e, deriv = TRUE)
    expect_equal(c(crossprod(X, m$dl_eta), sum(m$dl_extra)),
                 numDeriv::grad(f, c(b, e), what = "l"), tolerance = 1e-6, info = mg)
    expect_equal(c(crossprod(X, m$du_eta), sum(m$du_extra)),
                 numDeriv::grad(f, c(b, e), what = "u"), tolerance = 1e-6, info = mg)
  }
})

test_that("joint fit with Weibull, inverse Gaussian and Beta margins recovers parameters", {
  skip_on_cran()
  for (cop in c("gaussian", "student")) {
    d <- sim_new(2000, cop, seed = 3)
    fit <- tricopula(list(y1 ~ x, y2 ~ x, y3 ~ x), d, margins = new_margins, copula = cop)
    expect_equal(fit$convergence, 0L)
    expect_equal(names(coef(fit))[c(3, 6, 9)],
                 c("margin1_log_shape", "margin2_log_scale", "margin3_log_precision"))
    z <- (coef(fit) - true_new(cop)) / sqrt(diag(vcov(fit)))
    expect_lt(max(abs(z)), 4)
    expect_lt(max(abs(fit$gradient)), 1e-2)
    par <- true_new(cop) * 1.05
    expect_equal(unname(tc_gradient(fit, par)),
                 numDeriv::grad(function(p) tc_loglik(fit, p), par), tolerance = 1e-5)
  }
  expect_true(all(c("margin1_shape", "margin2_sigma", "margin3_precision") %in%
                    rownames(tc_dependence(fit))))
  pr <- predict(fit, data.frame(x = 0))
  expect_named(pr$margin1, c("scale", "shape"))
  expect_named(pr$margin3, c("mean", "precision"))
  s <- simulate(fit, seed = 4)
  expect_true(all(s$y1 > 0 & s$y2 > 0 & s$y3 > 0 & s$y3 < 1))
  mu <- predict(fit, type = "mean")
  expect_equal(colMeans(mu), colMeans(d[c("y1", "y2", "y3")]), tolerance = 0.05)
})

test_that("near-singular correlations are rejected instead of yielding a spurious likelihood", {
  skip_on_cran()
  # Regression test for run 6: with these data the copula step of the start
  # values wandered to |rho| = 1, where the density was rounding noise and
  # the reported log-likelihood was about +14900 instead of about -1357.
  d <- sim_new(1000, "student", seed = 18033)
  fit <- tricopula(list(y1 ~ x, y2 ~ x, y3 ~ x), d, margins = new_margins,
                   copula = "student")
  expect_equal(fit$convergence, 0L)
  expect_false(is.null(fit$vcov))
  expect_lt(fit$logLik, 0)
  expect_lt(max(abs(fit$correlation[upper.tri(diag(3))])), 0.9)
  expect_lt(max(abs(fit$gradient)), 1e-3)
  expect_lt(abs(fit$logLik - fit$ifm$logLik), 5)

  U <- matrix(runif(30), 10, 3)
  W <- rep(list(matrix(1, 10, 1, dimnames = list(NULL, "(Intercept)"))), 3)
  expect_true(all(TriCopula:::.tc_cop_eval(c(17.8, -17.6, 22.6), U, W, "gaussian")$lc == -Inf))
  expect_true(all(is.finite(TriCopula:::.tc_cop_eval(c(5.9, -5.9, 5.9), U, W, "gaussian")$lc)))
  expect_true(all(TriCopula:::.tc_pair(U[, 1], U[, 2], 20, "gaussian")$l == -Inf))
})

test_that("new margins accept scale covariates and reject values outside the support", {
  d <- sim_new(400, seed = 5)
  fit <- tricopula(list(y1 ~ x, y2 ~ x, y3 ~ x), d, margins = new_margins,
                   scale = ~ x, hessian = FALSE)
  expect_length(coef(fit), 15L)
  expect_true("margin3_log_precision_x" %in% names(coef(fit)))
  d$y3[1] <- 1.2
  expect_error(tricopula(list(y1 ~ x, y2 ~ x, y3 ~ x), d, margins = new_margins),
               "between 0 and 1")
  d$y3[1] <- 0.5; d$y1[2] <- -1
  expect_error(tricopula(list(y1 ~ x, y2 ~ x, y3 ~ x), d, margins = new_margins),
               "strictly positive")
})
