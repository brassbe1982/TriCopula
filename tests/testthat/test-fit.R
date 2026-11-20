test_that("analytic gradient matches numerical differentiation", {
  for (cop in c("gaussian", "student")) {
    d <- sim_mixed(200, cop, seed = 11)
    fit <- fit_mixed(d, cop, hessian = FALSE)
    set.seed(12)
    par <- true_par(cop) + rnorm(length(true_par(cop)), sd = 0.1)
    num <- numDeriv::grad(function(p) tc_loglik(fit, p), par)
    expect_equal(unname(tc_gradient(fit, par)), num, tolerance = 1e-5)
  }
})

test_that("per-observation scores sum to the gradient and are centred at the MLE", {
  d <- sim_mixed(300, seed = 13)
  fit <- fit_mixed(d)
  expect_equal(colSums(fit$scores), fit$gradient)
  expect_lt(max(abs(fit$gradient)), 1e-3)
})

test_that("fit estimates three continuous margins (0.1.0 smoke test)", {
  set.seed(101)
  d <- data.frame(x = rnorm(300))
  d$y1 <- 1 + .7 * d$x + rnorm(300)
  d$y2 <- exp(.3 + .2 * d$x + rnorm(300, sd = .4))
  d$y3 <- rgamma(300, shape = 3, scale = exp(.1 + .3 * d$x) / 3)
  fit <- tricopula(
    list(y1 ~ x, y2 ~ x, y3 ~ x), d,
    margins = c("gaussian", "lognormal", "gamma"),
    copula = "gaussian"
  )
  expect_s3_class(fit, "tricopula")
  expect_true(is.finite(as.numeric(logLik(fit))))
  expect_equal(dim(fit$correlation), c(3L, 3L))
  expect_gt(min(eigen(fit$correlation, symmetric = TRUE)$values), 0)
  expect_equal(length(tc_pit(fit)), 3L)
})

test_that("Gaussian margins + Gaussian copula reproduce the multivariate normal MLE", {
  skip_if_not_installed("mvtnorm")
  # With identical regressors in every equation the joint normal MLE is
  # equation-by-equation OLS with the ML residual covariance matrix.
  set.seed(21)
  n <- 400
  d <- data.frame(x = rnorm(n), w = runif(n))
  E <- matrix(rnorm(3 * n), n, 3) %*% chol(R_true) %*% diag(c(1, 2, 0.5))
  d$y1 <- 1 + 0.5 * d$x - d$w + E[, 1]
  d$y2 <- -1 + 0.2 * d$x + 2 * d$w + E[, 2]
  d$y3 <- 0.3 * d$x + E[, 3]
  fit <- tricopula(list(y1 ~ x + w, y2 ~ x + w, y3 ~ x + w), d,
                   margins = rep("gaussian", 3))
  ols <- lm(cbind(y1, y2, y3) ~ x + w, d)
  S <- crossprod(residuals(ols)) / n
  b <- coef(fit)
  expect_equal(unname(b[c(1:3, 5:7, 9:11)]), unname(c(coef(ols))), tolerance = 1e-5)
  expect_equal(unname(exp(b[c(4, 8, 12)])), unname(sqrt(diag(S))), tolerance = 1e-5)
  expect_equal(fit$correlation, unname(cov2cor(S)), tolerance = 1e-5)
  ll <- sum(mvtnorm::dmvnorm(residuals(ols), sigma = S, log = TRUE))
  expect_equal(as.numeric(logLik(fit)), ll, tolerance = 1e-8)
})

test_that("full ML recovers known parameters with valid standard errors", {
  skip_on_cran()
  for (cop in c("gaussian", "student")) {
    d <- sim_mixed(2000, cop, seed = 31)
    fit <- fit_mixed(d, cop)
    expect_equal(fit$convergence, 0L)
    expect_false(is.null(fit$vcov))
    se <- sqrt(diag(vcov(fit)))
    z <- (coef(fit) - true_par(cop)) / se
    expect_lt(max(abs(z)), 4)
    expect_equal(fit$correlation, R_true, tolerance = 0.1)
    expect_gte(fit$logLik, fit$ifm$logLik - 1e-6)
    dg <- tc_diagnostics(fit)
    expect_true(dg$hessian_positive_definite)
    expect_lt(dg$max_abs_gradient, 1e-2)
    expect_output(print(dg), "TriCopula diagnostics")
    expect_output(print(summary(fit)), "natural scale")
    expect_equal(dim(vcov(fit, type = "sandwich")), dim(vcov(fit)))
  }
})

test_that("Student-t fit to Gaussian-copula data warns and keeps usable standard errors", {
  hit <- FALSE
  for (seed in 200001:200008) {
    d <- sim_mixed(400, "gaussian", seed = seed)
    w <- NULL
    fit <- withCallingHandlers(
      fit_mixed(d, "student"),
      warning = function(cnd) { w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning") })
    if (fit$df > 100) {
      hit <- TRUE
      expect_match(paste(w, collapse = " "), "practically a Gaussian copula")
      expect_false(is.null(fit$vcov))
      se <- sqrt(diag(vcov(fit)))
      expect_true(all(is.finite(se[1:12])))
      if (fit$df_boundary) {
        expect_true(is.na(se[13]))
        expect_true(is.na(tc_dependence(fit)["df", "std.error"]))
        expect_true(all(is.finite(sqrt(diag(vcov(fit, type = "sandwich")))[1:12])))
      }
      expect_true(all(is.finite(tc_dependence(fit)[c("rho12", "rho13", "rho23"), "std.error"])))
      expect_output(print(summary(fit)), "natural scale")
      break
    }
  }
  expect_true(hit)
})

test_that("multiple starts agree and input errors are informative", {
  d <- sim_mixed(300, seed = 41)
  set.seed(1)
  f1 <- fit_mixed(d, hessian = FALSE)
  f3 <- fit_mixed(d, hessian = FALSE, nstart = 3)
  expect_equal(f1$logLik, f3$logLik, tolerance = 1e-6)
  expect_length(f3$starts$logLik, 3L)
  f_user <- fit_mixed(d, hessian = FALSE, start = true_par())
  expect_equal(f_user$logLik, f1$logLik, tolerance = 1e-6)
  expect_null(f_user$ifm)
  expect_error(fit_mixed(d, start = 1:3), "incorrect length")
  d$y2[1] <- -1
  expect_error(fit_mixed(d), "strictly positive")
  d$y2[1] <- NA
  expect_error(fit_mixed(d))
  expect_error(tricopula(list(y1 ~ x, y2 ~ x), d), "three formulas")
})
