delta_true <- list(c(0.5, 0.3), c(0.3, -0.2), c(0.3, 0.2))

sim_dep <- function(n, copula = "gaussian", df = 5, seed = 1) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n))
  U <- TriCopula:::.tc_rcopula_vec(
    n, tanh(delta_true[[1]][1] + delta_true[[1]][2] * d$x),
    tanh(delta_true[[2]][1] + delta_true[[2]][2] * d$x),
    tanh(delta_true[[3]][1] + delta_true[[3]][2] * d$x), copula, df)
  d$y1 <- qnorm(U[, 1], 1 + 0.5 * d$x, 1.2)
  d$y2 <- qlnorm(U[, 2], 0.2 + 0.3 * d$x, 0.5)
  d$y3 <- qgamma(U[, 3], shape = 2.5, scale = exp(0.1 + 0.2 * d$x) / 2.5)
  d
}

true_dep_par <- function(copula = "gaussian", df = 5) {
  th <- c(1, 0.5, log(1.2), 0.2, 0.3, log(0.5), 0.1, 0.2, log(2.5),
          unlist(delta_true))
  if (copula == "student") th <- c(th, log(df - 2.05))
  th
}

test_that("observation-specific correlation matrices are valid and the vectorized density is correct", {
  skip_if_not_installed("mvtnorm")
  set.seed(1)
  n <- 25
  dep <- TriCopula:::.tc_dep_vec(rnorm(n, sd = 2), rnorm(n, sd = 2), rnorm(n, sd = 2))
  U <- matrix(runif(3 * n), n, 3)
  for (cop in c("gaussian", "student")) {
    df <- if (cop == "student") 4.5 else Inf
    lc <- TriCopula:::.tc_logcop(U, dep, df, cop)$lc
    ref <- vapply(seq_len(n), function(i) {
      R <- TriCopula:::.tc_R(dep$r12[i], dep$r13[i], dep$r23[i])
      expect_gt(min(eigen(R, symmetric = TRUE)$values), 0)
      TriCopula:::.tc_logcop(U[i, , drop = FALSE], R, df, cop)$lc
    }, numeric(1))
    expect_equal(lc, ref, tolerance = 1e-9)
    z <- if (cop == "gaussian") qnorm(U) else qt(U, df)
    ref2 <- vapply(seq_len(n), function(i) {
      R <- TriCopula:::.tc_R(dep$r12[i], dep$r13[i], dep$r23[i])
      if (cop == "gaussian") {
        mvtnorm::dmvnorm(z[i, ], sigma = R, log = TRUE) - sum(dnorm(z[i, ], log = TRUE))
      } else {
        mvtnorm::dmvt(z[i, ], sigma = R, df = df, log = TRUE) - sum(dt(z[i, ], df, log = TRUE))
      }
    }, numeric(1))
    expect_equal(lc, ref2, tolerance = 1e-8)
  }
})

test_that("the vectorized sampler reproduces the target correlations", {
  set.seed(2)
  n <- 2e5
  U <- TriCopula:::.tc_rcopula_vec(n, rep(0.5, n), rep(0.3, n),
                                   rep((0.4 - 0.15) / sqrt(0.75 * 0.91), n),
                                   "gaussian", Inf)
  expect_equal(cor(qnorm(U)), R_true, tolerance = 0.01, ignore_attr = TRUE)
})

test_that("intercept-only dependence reproduces the constant-dependence fit", {
  d <- sim_mixed(300, seed = 71)
  f0 <- fit_mixed(d)
  f1 <- fit_mixed(d, dependence = list(~ 1, ~ 1, ~ 1))
  expect_equal(coef(f0), coef(f1))
  expect_true(f0$dep_constant)
  expect_equal(names(coef(f0))[10:12],
               c("atanh_partial_rho12", "atanh_partial_rho13",
                 "atanh_partial_rho23_given_1"))
  expect_equal(nrow(f0$rho), 300L)
  expect_equal(unname(f0$rho[1, 1:3]), f0$correlation[c(4, 7, 8)])
})

test_that("gradient is correct with dependence covariates", {
  for (cop in c("gaussian", "student")) {
    d <- sim_dep(200, cop, seed = 72)
    fit <- fit_mixed(d, cop, hessian = FALSE, dependence = ~ x)
    set.seed(73)
    par <- true_dep_par(cop) + rnorm(length(true_dep_par(cop)), sd = 0.1)
    num <- numDeriv::grad(function(p) tc_loglik(fit, p), par)
    expect_equal(unname(tc_gradient(fit, par)), num, tolerance = 1e-5)
  }
})

test_that("dependence coefficients are recovered and post-estimation follows them", {
  skip_on_cran()
  d <- sim_dep(3000, seed = 74)
  fit <- fit_mixed(d, dependence = ~ x)
  expect_equal(fit$convergence, 0L)
  expect_false(fit$dep_constant)
  expect_length(coef(fit), 15L)
  expect_true("atanh_partial_rho12_x" %in% names(coef(fit)))
  z <- (coef(fit) - true_dep_par()) / sqrt(diag(vcov(fit)))
  expect_lt(max(abs(z)), 4)
  expect_gte(fit$logLik, fit$ifm$logLik - 1e-6)

  f0 <- fit_mixed(d)
  lr <- tc_lrtest(f0, fit)
  expect_equal(lr$df, 3)
  expect_lt(lr$p.value, 1e-6)
  expect_error(tc_lrtest(fit, f0), "more parameters")

  nd <- data.frame(x = c(-1, 0, 1))
  pd <- predict(fit, nd, type = "dependence")
  b <- coef(fit)
  expect_equal(pd$rho12, unname(tanh(b["atanh_partial_rho12_(Intercept)"] +
                                       b["atanh_partial_rho12_x"] * nd$x)))
  expect_equal(pd$rho23, pd$rho12 * pd$rho13 +
                 sqrt((1 - pd$rho12^2) * (1 - pd$rho13^2)) * pd[["rho23|1"]])
  expect_equal(pd$tau12, 2 / pi * asin(pd$rho12))
  expect_equal(dim(predict(fit, type = "dependence")), c(3000L, 7L))

  # average-correlation summaries
  dep <- tc_dependence(fit)
  expect_equal(dep["rho12", "estimate"], mean(fit$rho[, "rho12"]))
  expect_equal(fit$correlation[1, 2], mean(fit$rho[, "rho12"]))
  expect_true(all(dep$std.error > 0))

  # simulation honours row-specific dependence
  set.seed(75)
  big <- data.frame(x = rep(c(-1.5, 1.5), each = 4e4))
  s <- simulate(fit, newdata = big)
  r_lo <- cor(s$y1[big$x < 0], log(s$y2[big$x < 0]))
  r_hi <- cor(s$y1[big$x > 0], log(s$y2[big$x > 0]))
  pr <- predict(fit, data.frame(x = c(-1.5, 1.5)), type = "dependence")$rho12
  expect_equal(c(r_lo, r_hi), pr, tolerance = 0.03)

  # joint probability honours row-specific dependence
  skip_if_not_installed("mvtnorm")
  th <- c(1, 1.2, 1)
  jp <- tc_jointprob(fit, th, data.frame(x = 1.5))
  mc <- with(s[big$x > 0, ], mean(y1 <= th[1] & y2 <= th[2] & y3 <= th[3]))
  expect_equal(jp, mc, tolerance = 0.02)
})

test_that("different formulas per dependence parameter and input errors", {
  d <- sim_dep(400, "student", seed = 76)
  d$g <- factor(rep(c("a", "b"), 200))
  fit <- fit_mixed(d, "student", dependence = list(~ x, ~ 1, ~ g))
  expect_equal(fit$index$q, c(2L, 1L, 2L))
  expect_true(all(c("atanh_partial_rho13", "atanh_partial_rho23_given_1_gb")
                  %in% names(coef(fit))))
  pd <- predict(fit, data.frame(x = 0, g = factor("b", levels = c("a", "b"))),
                type = "dependence")
  expect_equal(nrow(pd), 1L)
  expect_error(fit_mixed(d, dependence = y1 ~ x), "one-sided")
  expect_error(fit_mixed(d, dependence = list(~ x, ~ 1)), "list of three")
})
