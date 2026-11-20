# Covariate-dependent degrees of freedom (0.10.0)

sim_dft <- function(n, seed) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), g = rbinom(n, 1, 0.5))
  dfi <- 2.05 + exp(log(2) + log(6) * d$g)        # df 4.05 (g = 0) or 14.05 (g = 1)
  Z <- matrix(rnorm(3 * n), n, 3) %*% chol(R_true)
  U <- pt(Z / sqrt(rchisq(n, dfi) / dfi), dfi)
  d$y1 <- qnorm(U[, 1], 1 + 0.5 * d$x, 1.2)
  d$y2 <- qlnorm(U[, 2], 0.2 + 0.3 * d$x, 0.5)
  d$y3 <- qgamma(U[, 3], shape = 2.5, scale = exp(0.1 + 0.2 * d$x) / 2.5)
  d
}

test_that("df = ~ 1 is the default model and a zero slope is the constant-df model", {
  skip_on_cran()
  d <- sim_mixed(300, "student", seed = 141)
  f0 <- fit_mixed(d, "student", hessian = FALSE)
  f1 <- fit_mixed(d, "student", hessian = FALSE, df = ~ 1)
  expect_identical(coef(f1), coef(f0))
  expect_identical(f1$logLik, f0$logLik)
  expect_true(f0$df_constant)
  expect_null(f0$df_i); expect_null(f0$df_design)
  expect_equal(names(coef(f0))[13], "log_df_minus_2.05")
  f2 <- fit_mixed(d, "student", hessian = FALSE, df = ~ x, start = c(coef(f0), 0),
                  control = list(maxit = 0), polish = FALSE)
  expect_equal(names(coef(f2))[13:14],
               c("log_df_minus_2.05_(Intercept)", "log_df_minus_2.05_x"))
  expect_equal(tc_loglik(f2, c(coef(f0), 0)), f0$logLik, tolerance = 1e-10)
  expect_equal(unname(tc_gradient(f2, c(coef(f0), 0))[1:13]), unname(f0$gradient),
               tolerance = 1e-8)
  # PCC: same for a Student-t pair
  fam <- c("student", "clayton", "gaussian")
  p0 <- fit_mixed(d, "pcc", families = fam, hessian = FALSE)
  p1 <- fit_mixed(d, "pcc", families = fam, hessian = FALSE, df = ~ 1, root = 1)
  expect_identical(coef(p1), coef(p0))
  p2 <- fit_mixed(d, "pcc", families = fam, hessian = FALSE, df = ~ x,
                  start = c(coef(p0), 0), control = list(maxit = 0), polish = FALSE)
  expect_equal(names(coef(p2))[13:14],
               c("pcc12_log_df_minus_2.05_(Intercept)", "pcc12_log_df_minus_2.05_x"))
  expect_equal(tc_loglik(p2, c(coef(p0), 0)), p0$logLik, tolerance = 1e-10)
  expect_error(fit_mixed(d, "gaussian", df = ~ x), "only used")
  expect_error(fit_mixed(d, "pcc", families = c("clayton", "gumbel", "frank"), df = ~ x),
               "only used")
  expect_error(fit_mixed(d, "student", df = list(~ x, ~ 1, ~ 1)), "single")
})

test_that("log-likelihood with covariate-dependent df matches mvtnorm and VineCopula", {
  skip_if_not_installed("mvtnorm")
  d <- sim_dft(120, 142)
  par <- c(true_par("student"), 0.4)
  par[13] <- log(3)
  f <- fit_mixed(d, "student", hessian = FALSE, df = ~ x, start = par,
                 control = list(maxit = 0), polish = FALSE)
  dfi <- 2.05 + exp(par[13] + par[14] * d$x)
  u <- cbind(pnorm(d$y1, 1 + 0.5 * d$x, 1.2), plnorm(d$y2, 0.2 + 0.3 * d$x, 0.5),
             pgamma(d$y3, shape = 2.5, scale = exp(0.1 + 0.2 * d$x) / 2.5))
  lm <- dnorm(d$y1, 1 + 0.5 * d$x, 1.2, log = TRUE) +
    dlnorm(d$y2, 0.2 + 0.3 * d$x, 0.5, log = TRUE) +
    dgamma(d$y3, shape = 2.5, scale = exp(0.1 + 0.2 * d$x) / 2.5, log = TRUE)
  lc <- vapply(seq_len(nrow(d)), function(i) {
    z <- qt(u[i, ], dfi[i])
    mvtnorm::dmvt(z, sigma = R_true, df = dfi[i], log = TRUE) - sum(dt(z, dfi[i], log = TRUE))
  }, numeric(1))
  expect_equal(tc_loglik(f, par), sum(lm + lc), tolerance = 1e-9)
  expect_equal(unname(tc_joint_logscore(f)), lm + lc, tolerance = 1e-9)
  expect_equal(unname(tc_gradient(f, par)),
               numDeriv::grad(function(p) tc_loglik(f, p), par), tolerance = 1e-5)
  expect_equal(f$df_i, dfi, tolerance = 1e-12)
  expect_equal(f$df, mean(dfi), tolerance = 1e-12)

  skip_if_not_installed("VineCopula")
  fam <- c("student", "joe", "student")
  pp <- c(true_par()[1:9], atanh(0.5), log(0.7), atanh(0.3),
          log(4), 0.5, log(6), -0.3)
  g <- fit_mixed(d, "pcc", families = fam, hessian = FALSE, df = list(~ x, ~ 1, ~ x),
                 start = pp, control = list(maxit = 0), polish = FALSE)
  expect_equal(names(coef(g))[13:16],
               c("pcc12_log_df_minus_2.05_(Intercept)", "pcc12_log_df_minus_2.05_x",
                 "pcc23_given_1_log_df_minus_2.05_(Intercept)",
                 "pcc23_given_1_log_df_minus_2.05_x"))
  d1 <- 2.05 + exp(log(4) + 0.5 * d$x); d3 <- 2.05 + exp(log(6) - 0.3 * d$x)
  n <- nrow(d)
  l12 <- log(VineCopula::BiCopPDF(u[, 1], u[, 2], rep(2, n), rep(0.5, n), d1))
  l13 <- log(VineCopula::BiCopPDF(u[, 1], u[, 3], 6, 1.7))
  a <- VineCopula::BiCopHfunc1(u[, 1], u[, 2], rep(2, n), rep(0.5, n), d1)
  b <- VineCopula::BiCopHfunc1(u[, 1], u[, 3], 6, 1.7)
  l23 <- log(VineCopula::BiCopPDF(a, b, rep(2, n), rep(0.3, n), d3))
  expect_equal(tc_loglik(g, pp), sum(lm + l12 + l13 + l23), tolerance = 1e-8)
  expect_equal(unname(tc_gradient(g, pp)),
               numDeriv::grad(function(p) tc_loglik(g, p), pp), tolerance = 1e-4)
  expect_equal(g$df_i[, 1], d1, tolerance = 1e-12)
  expect_true(all(is.na(g$df_i[, 2])))
})

test_that("fit with covariate-dependent df: recovery and post-estimation", {
  skip_on_cran()
  d <- sim_dft(3000, 143)
  fit <- fit_mixed(d, "student", df = ~ g)
  expect_equal(fit$convergence, 0L)
  expect_false(is.null(fit$vcov))
  truth <- c(true_par(), log(2), log(6))
  z <- (coef(fit) - truth) / sqrt(diag(vcov(fit)))
  expect_lt(max(abs(z)), 4)
  expect_lt(max(abs(fit$gradient)), 1e-2)
  expect_false(fit$df_constant)
  expect_gt(coef(fit)[["log_df_minus_2.05_g"]], 0)
  f0 <- fit_mixed(d, "student")
  lr <- tc_lrtest(f0, fit)
  expect_equal(lr$df, 1)
  expect_lt(lr$p.value, 0.01)
  pd <- predict(fit, data.frame(x = 0, g = 0:1), type = "dependence")
  expect_true("df" %in% names(pd))
  expect_equal(pd$df, 2.05 + exp(coef(fit)[[13]] + coef(fit)[[14]] * 0:1), tolerance = 1e-12)
  expect_false("df" %in% names(predict(f0, data.frame(x = 0, g = 0), type = "dependence")))
  dep <- tc_dependence(fit)
  expect_equal(dep["df", "estimate"], mean(fit$df_i), tolerance = 1e-10)
  expect_true(is.finite(dep["df", "std.error"]))
  expect_output(print(fit), "Average copula degrees of freedom")
  expect_s3_class(tc_diagnostics(fit), "tc_diagnostics")
  expect_true(all(is.finite(tc_taildep(fit))))
  bt <- tc_bootstrap(fit, B = 4, seed = 3)
  expect_equal(bt$failures, 0L)
  expect_equal(nrow(bt$coefficients), 14L)
  ef <- tc_effects(fit, "x")
  expect_equal(ef$effect[1], coef(fit)[["margin1_x"]], tolerance = 1e-6)
  skip_if_not_installed("mvtnorm")
  # joint probabilities differ by group and match simulation from the fit
  th <- c(-0.5, 0.6, 0.4)
  jp <- tc_jointprob(fit, th, data.frame(x = 0, g = 0:1))
  set.seed(144)
  for (gv in 0:1) {
    s <- simulate(fit, newdata = data.frame(x = rep(0, 1e5), g = gv))
    expect_equal(jp[gv + 1], mean(s$y1 <= th[1] & s$y2 <= th[2] & s$y3 <= th[3]),
                 tolerance = 0.03)
  }
  expect_gt(jp[1], jp[2])      # heavier joint tail with the smaller df
})

test_that("PCC fit with a covariate-dependent df of a Student-t pair", {
  skip_on_cran()
  set.seed(145)
  n <- 2500
  d <- data.frame(x = rnorm(n), g = rbinom(n, 1, 0.5))
  U <- TriCopula:::.tc_rpcc(n, list(atanh(0.6), log(1.5), atanh(0.3)),
                            c("student", "clayton", "gaussian"),
                            list(log(1.5) + log(8) * d$g, NULL, NULL))
  d$y1 <- qnorm(U[, 1], 1 + 0.5 * d$x, 1.2)
  d$y2 <- qlnorm(U[, 2], 0.2 + 0.3 * d$x, 0.5)
  d$y3 <- qgamma(U[, 3], shape = 2.5, scale = exp(0.1 + 0.2 * d$x) / 2.5)
  fam <- c("student", "clayton", "gaussian")
  fit <- fit_mixed(d, "pcc", families = fam, df = ~ g)
  expect_equal(fit$convergence, 0L)
  truth <- c(true_par()[1:9], atanh(0.6), log(1.5), atanh(0.3), log(1.5), log(8))
  z <- (coef(fit) - truth) / sqrt(diag(vcov(fit)))
  expect_lt(max(abs(z)), 4)
  expect_equal(names(coef(fit))[13:14],
               c("pcc12_log_df_minus_2.05_(Intercept)", "pcc12_log_df_minus_2.05_g"))
  pd <- predict(fit, data.frame(x = 0, g = 0:1), type = "dependence")
  expect_true("df12" %in% names(pd))
  expect_lt(pd$df12[1], pd$df12[2])
  expect_true("df12" %in% rownames(tc_dependence(fit)))
  s <- simulate(fit, seed = 1)
  expect_equal(dim(s), c(n, 3L))
  skip_if_not_installed("mvtnorm")
  jp <- tc_jointprob(fit, c(1, 1.2, 1), data.frame(x = 0, g = 0:1))
  expect_true(all(jp > 0 & jp < 1))
  expect_false(isTRUE(all.equal(jp[1], jp[2])))
})
