gamma_true <- list(c(log(1.2), 0.3), c(log(0.5), -0.2), c(log(2.5), 0.3))

sim_ls <- function(n, copula = "gaussian", df = 5, seed = 1, dep_x = FALSE) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n))
  U <- if (dep_x) {
    TriCopula:::.tc_rcopula_vec(n, tanh(0.5 + 0.3 * d$x), tanh(0.3 - 0.2 * d$x),
                                tanh(0.3 + 0.2 * d$x), copula, df)
  } else tc_rcopula(n, R_true, copula, df)
  d$y1 <- qnorm(U[, 1], 1 + 0.5 * d$x, exp(gamma_true[[1]][1] + gamma_true[[1]][2] * d$x))
  d$y2 <- qlnorm(U[, 2], 0.2 + 0.3 * d$x, exp(gamma_true[[2]][1] + gamma_true[[2]][2] * d$x))
  sh <- exp(gamma_true[[3]][1] + gamma_true[[3]][2] * d$x)
  d$y3 <- qgamma(U[, 3], shape = sh, scale = exp(0.1 + 0.2 * d$x) / sh)
  d
}

true_ls_par <- function(copula = "gaussian", df = 5, dep_x = FALSE) {
  th <- c(1, 0.5, gamma_true[[1]], 0.2, 0.3, gamma_true[[2]], 0.1, 0.2, gamma_true[[3]],
          if (dep_x) c(0.5, 0.3, 0.3, -0.2, 0.3, 0.2) else TriCopula:::.tc_dep_inverse(R_true))
  if (copula == "student") th <- c(th, log(df - 2.05))
  th
}

test_that("intercept-only scale reproduces the constant-scale fit and names", {
  d <- sim_mixed(300, seed = 81)
  f0 <- fit_mixed(d)
  f1 <- fit_mixed(d, scale = list(~ 1, ~ 1, ~ 1))
  expect_equal(coef(f0), coef(f1))
  expect_true(f0$scale_constant)
  expect_equal(names(coef(f0))[c(3, 6, 9)],
               c("margin1_log_scale", "margin2_log_scale", "margin3_log_shape"))
})

test_that("gradient is correct with scale and dependence covariates", {
  for (cop in c("gaussian", "student")) {
    d <- sim_ls(200, cop, seed = 82, dep_x = TRUE)
    fit <- fit_mixed(d, cop, hessian = FALSE, scale = ~ x, dependence = ~ x)
    set.seed(83)
    par <- true_ls_par(cop, dep_x = TRUE)
    par <- par + rnorm(length(par), sd = 0.1)
    num <- numDeriv::grad(function(p) tc_loglik(fit, p), par)
    expect_equal(unname(tc_gradient(fit, par)), num, tolerance = 1e-5)
  }
})

test_that("heteroscedastic Gaussian margin matches a direct normal likelihood", {
  # Independent responses: the joint fit with rho fixed near its estimate
  # must agree with a stand-alone location-scale normal MLE for margin 1.
  d <- sim_ls(500, seed = 84)
  X <- cbind(1, d$x)
  nll <- function(th) -sum(dnorm(d$y1, X %*% th[1:2], exp(X %*% th[3:4]), log = TRUE))
  ref <- optim(c(1, 0.5, 0, 0), nll, method = "BFGS", control = list(reltol = 1e-12))
  mj <- TriCopula:::.tc_fit_margin(d$y1, X, "gaussian",
                                   matrix(X, ncol = 2, dimnames = list(NULL, c("(Intercept)", "x"))))
  expect_equal(unname(mj$par), ref$par, tolerance = 1e-4)
  expect_equal(mj$logLik, -ref$value, tolerance = 1e-6)
})

test_that("scale coefficients are recovered and post-estimation follows them", {
  skip_on_cran()
  d <- sim_ls(3000, seed = 85)
  fit <- fit_mixed(d, scale = ~ x)
  expect_equal(fit$convergence, 0L)
  expect_false(fit$scale_constant)
  expect_length(coef(fit), 15L)
  expect_true(all(c("margin1_log_scale_x", "margin3_log_shape_x") %in% names(coef(fit))))
  z <- (coef(fit) - true_ls_par()) / sqrt(diag(vcov(fit)))
  expect_lt(max(abs(z)), 4)
  expect_gte(fit$logLik, fit$ifm$logLik - 1e-6)
  expect_lt(tc_lrtest(fit_mixed(d), fit)$p.value, 1e-6)

  nd <- data.frame(x = c(-1, 0, 1))
  b <- coef(fit)
  pr <- predict(fit, nd)
  expect_equal(pr$margin1$scale,
               unname(exp(b["margin1_log_scale_(Intercept)"] + b["margin1_log_scale_x"] * nd$x)))
  expect_equal(pr$margin3$shape,
               unname(exp(b["margin3_log_shape_(Intercept)"] + b["margin3_log_shape_x"] * nd$x)))
  mu <- predict(fit, nd, type = "mean")
  expect_equal(mu$y2, exp(pr$margin2$meanlog + pr$margin2$sdlog^2 / 2))

  dep <- tc_dependence(fit)
  expect_equal(dep["margin1_sigma", "estimate"],
               mean(exp(drop(fit$scale_design[[1]] %*% b[3:4]))))

  set.seed(86)
  big <- data.frame(x = rep(c(-1, 1), each = 5e4))
  s <- simulate(fit, newdata = big)
  expect_equal(c(sd(s$y1[big$x < 0]), sd(s$y1[big$x > 0])),
               predict(fit, data.frame(x = c(-1, 1)))$margin1$scale, tolerance = 0.02)

  skip_if_not_installed("mvtnorm")
  th <- c(1, 1.2, 1)
  jp <- tc_jointprob(fit, th, data.frame(x = 1))
  mc <- with(s[big$x > 0, ], mean(y1 <= th[1] & y2 <= th[2] & y3 <= th[3]))
  expect_equal(jp, mc, tolerance = 0.02)
})

test_that("per-margin scale formulas and input errors", {
  d <- sim_ls(400, seed = 87)
  fit <- fit_mixed(d, scale = list(~ x, ~ 1, ~ x), dependence = list(~ 1, ~ x, ~ 1))
  expect_equal(fit$index$e, c(2L, 1L, 2L))
  expect_equal(fit$index$q, c(1L, 2L, 1L))
  expect_true("margin2_log_scale" %in% names(coef(fit)))
  expect_lt(max(abs(fit$gradient)), 1e-2)
  expect_error(fit_mixed(d, scale = y1 ~ x), "one-sided")
  expect_error(fit_mixed(d, scale = list(~ x)), "list of three")
})
