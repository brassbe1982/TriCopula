# Non-simplified vine (0.15.0): the conditional pair copula depends on the
# conditioning value through a penalized spline

sim_ns <- function(n = 600, f = function(u) 1.2 * (u - 0.5), seed = 301,
                   fam = c("clayton", "gumbel", "gaussian")) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n))
  w <- matrix(runif(3 * n), n, 3)
  hinv <- TriCopula:::.tc_pair_hinv
  u1 <- w[, 1]
  u2 <- hinv(w[, 2], u1, log(1.2), fam[1])
  b <- hinv(w[, 3], w[, 2], 0.4 + f(u1), fam[3])
  u3 <- hinv(b, u1, log(0.5), fam[2])
  d$y1 <- qnorm(u1, 1 + 0.5 * d$x, 1.2)
  d$y2 <- qlnorm(u2, 0.2 + 0.3 * d$x, 0.5)
  d$y3 <- qgamma(u3, shape = 2.5, scale = exp(0.1 + 0.2 * d$x) / 2.5)
  d
}
fam_ns <- c("clayton", "gumbel", "gaussian")
fm_ns <- list(y1 ~ x, y2 ~ x, y3 ~ x)

test_that("defaults unchanged: simplified = TRUE is the earlier model", {
  d <- sim_ns(300)
  f <- tricopula(fm_ns, d, copula = "pcc", families = fam_ns, hessian = FALSE)
  g <- tricopula(fm_ns, d, copula = "pcc", families = fam_ns, hessian = FALSE,
                 simplified = TRUE, k_cond = 9)
  expect_identical(coef(f), coef(g))
  expect_identical(f$logLik, g$logLik)
  expect_null(f$nonsimplified)
  expect_true(f$simplified)
  expect_length(f$index$ns, 0L)
  expect_error(tricopula(fm_ns, d, copula = "gaussian", simplified = FALSE), "only used with")
  expect_error(tricopula(fm_ns, d, copula = "pcc", families = fam_ns, simplified = NA),
               "TRUE or FALSE")
  expect_error(tricopula(fm_ns, d, copula = "pcc", families = fam_ns, simplified = FALSE,
                         k_cond = 3), "k_cond")
})

test_that("spline basis: constraint, derivative and penalty", {
  ns <- TriCopula:::.tc_ns_setup(7)
  expect_equal(ns$ncoef, 6L)
  # every basis function integrates to zero over (0, 1)
  for (j in 1:6) {
    expect_lt(abs(integrate(function(u) ns$basis(u, FALSE)$B[, j], 0, 1,
                            rel.tol = 1e-10)$value), 1e-8)
  }
  u <- c(0.03, 0.2, 0.5, 0.77, 0.98)
  g <- c(0.5, -1, 0.3, 0.8, -0.2, 0.4)
  fd <- numDeriv::grad(function(v) sum(ns$basis(v, FALSE)$B %*% g), u)
  expect_equal(drop(ns$basis(u)$dB %*% g), fd, tolerance = 1e-7)
  # second-difference penalty: positive semi-definite, linear functions unpenalized
  ev <- eigen(ns$S, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(ev > -1e-10))
  expect_equal(sum(ev > 1e-8), 5L)
  # the constant is not in the constrained basis, the centred linear function is
  uu <- seq(0.01, 0.99, length.out = 200)
  B <- ns$basis(uu, FALSE)$B
  expect_lt(max(abs(lm.fit(B, uu - 0.5)$residuals)), 1e-8)
  expect_gt(max(abs(lm.fit(B, rep(1, 200))$residuals)), 0.1)
})

test_that("log-density of the non-simplified vine equals VineCopula with per-row parameters", {
  skip_if_not_installed("VineCopula")
  set.seed(5)
  n <- 200
  U <- matrix(runif(3 * n, 0.02, 0.98), n, 3)
  ns <- TriCopula:::.tc_ns_setup(6)
  W <- rep(list(matrix(1, n, 1, dimnames = list(NULL, "(Intercept)"))), 3)
  g <- c(0.6, -0.4, 0.5, -0.7, 0.3)
  theta <- c(log(1.3), log(0.6), 0.3, g)
  ev <- TriCopula:::.tc_pcc_eval(theta, U, W, fam_ns, ns = ns)
  eta3 <- 0.3 + drop(ns$basis(U[, 1], FALSE)$B %*% g)
  h2 <- VineCopula::BiCopHfunc1(U[, 1], U[, 2], family = 3, par = 1.3)
  h3 <- VineCopula::BiCopHfunc1(U[, 1], U[, 3], family = 4, par = 1.6)
  ref <- log(VineCopula::BiCopPDF(U[, 1], U[, 2], family = 3, par = 1.3)) +
    log(VineCopula::BiCopPDF(U[, 1], U[, 3], family = 4, par = 1.6)) +
    log(VineCopula::BiCopPDF(h2, h3, family = rep(1, n), par = tanh(eta3)))
  expect_equal(ev$lc, ref, tolerance = 1e-8)
  expect_equal(ev$dep$eta[[3]], eta3)
  # with zero spline coefficients it is the simplified vine
  ev0 <- TriCopula:::.tc_pcc_eval(c(theta[1:3], rep(0, 5)), U, W, fam_ns, ns = ns)
  expect_equal(ev0$lc, TriCopula:::.tc_pcc_eval(theta[1:3], U, W, fam_ns)$lc)
})

test_that("non-simplified fit: score, nesting, test, smooth and post-estimation", {
  skip_on_cran()
  d <- sim_ns(600)
  f0 <- tricopula(fm_ns, d, copula = "pcc", families = fam_ns)
  f1 <- tricopula(fm_ns, d, copula = "pcc", families = fam_ns, simplified = FALSE,
                  k_cond = 6)
  expect_equal(f1$convergence, 0L)
  expect_equal(names(coef(f1))[1:12], names(coef(f0)))
  expect_equal(names(coef(f1))[13:17], paste0("pcc23_given_1_s(u1).", 1:5))
  # analytic score against numerical derivatives away from the optimum
  p <- coef(f1) + 0.05 * sin(seq_along(coef(f1)))
  expect_equal(unname(tc_gradient(f1, p)),
               numDeriv::grad(function(q) tc_loglik(f1, q), p), tolerance = 1e-5)
  expect_lt(max(abs(f1$gradient)), 1e-3)
  # the simplified model is the special case with zero spline coefficients
  expect_equal(tc_loglik(f1, c(coef(f0), rep(0, 5))), f0$logLik, tolerance = 1e-8)
  expect_gt(f1$logLik, f0$logLik)
  expect_lt(AIC(f1), AIC(f0))
  # test of the simplifying assumption (truth: linear in u1, slope 1.2)
  tt <- tc_simplifying_test(f0)
  expect_equal(tt$df, 3L)
  expect_equal(tt$statistic, 2 * (attr(tt, "fit")$logLik - f0$logLik))
  expect_lt(tt$p.value, 0.001)
  expect_error(tc_simplifying_test(f1), "simplified fit")
  # the unpenalized alternative fits at least as well as the penalized one of the same size
  f1u <- tricopula(fm_ns, d, copula = "pcc", families = fam_ns, simplified = FALSE,
                   k_cond = 4, sp = 0, hessian = FALSE)
  expect_equal(attr(tt, "fit")$logLik, f1u$logLik, tolerance = 1e-4)
  # fitted curve
  st <- summary(f1)$smooths
  expect_equal(st$term, "dep23|1:s(u1)")
  expect_true(st$edf > 0.9 && st$edf <= 5)
  sm <- tc_smooth(f1, "dep23|1:s(u1)", newdata = data.frame(u = c(0.1, 0.5, 0.9)))
  expect_equal(sm$fit, drop(f1$nonsimplified$basis(c(0.1, 0.5, 0.9))$B %*% coef(f1)[13:17]))
  expect_gt(sm$fit[3] - sm$fit[1], 0.4)                 # increasing, as the truth
  expect_true(all(sm$lower < sm$fit & sm$fit < sm$upper))
  expect_error(tc_smooth(f1, 1, newdata = data.frame(u = 2)), "values in")
  # per-observation dependence and its averages
  u1 <- f1$uniforms[, 1]
  eta3 <- coef(f1)[[12]] + drop(f1$nonsimplified$basis(u1)$B %*% coef(f1)[13:17])
  expect_equal(unname(f1$pair[, 3]), tanh(eta3), tolerance = 1e-10)
  dn <- tc_dependence(f1)
  expect_equal(dn["par23|1", "estimate"], mean(tanh(eta3)), tolerance = 1e-10)
  expect_equal(predict(f1, type = "dependence")[["par23|1"]], unname(f1$pair[, 3]))
  # bootstrap refits keep the model
  bt <- tc_bootstrap(f1, B = 2, seed = 1)
  expect_equal(colnames(bt$replicates)[1:17], names(coef(f1)))
})

test_that("non-simplified vine: sampler, joint probability and missing root", {
  skip_on_cran()
  d <- sim_ns(400, f = function(u) 0.8 * cos(2 * pi * u), seed = 302)
  f1 <- tricopula(fm_ns, d, copula = "pcc", families = fam_ns, simplified = FALSE,
                  k_cond = 6, hessian = FALSE)
  # joint probability against simulation from the fitted model
  thr <- c(1, 1.2, 1)
  set.seed(9)
  s <- simulate(f1, newdata = data.frame(x = rep(0, 200000)))
  emp <- mean(s[[1]] <= thr[1] & s[[2]] <= thr[2] & s[[3]] <= thr[3])
  expect_equal(tc_jointprob(f1, thr, data.frame(x = 0)), emp, tolerance = 0.02)
  # the simulated dependence of the conditional pair varies with u1 as fitted
  u <- sapply(1:3, function(j) rank(s[[j]]) / (nrow(s) + 1))
  h2 <- TriCopula:::.tc_pair(u[, 1], u[, 2], coef(f1)[[10]], fam_ns[1])$h
  h3 <- TriCopula:::.tc_pair(u[, 1], u[, 3], coef(f1)[[11]], fam_ns[2])$h
  mid <- u[, 1] > 0.4 & u[, 1] < 0.6
  lowhigh <- u[, 1] < 0.1 | u[, 1] > 0.9
  expect_gt(cor(qnorm(h2[lowhigh]), qnorm(h3[lowhigh])) - cor(qnorm(h2[mid]), qnorm(h3[mid])), 0.5)
  # missing root: the integral over the root against stats::integrate
  ns <- f1$nonsimplified
  cop <- TriCopula:::.tc_cop(f1)
  th <- coef(f1)[10:17]
  W1 <- lapply(f1$dep_design, function(w) w[1, , drop = FALSE])
  Um <- matrix(c(NA, 0.3, 0.8), 1)
  val <- TriCopula:::.tc_cop_eval_missing(th, Um, W1, cop)$lc
  ref <- integrate(function(s) {
    vapply(s, function(si) {
      exp(TriCopula:::.tc_pcc_eval(th, matrix(c(si, 0.3, 0.8), 1), W1, fam_ns, ns = ns)$lc)
    }, numeric(1))
  }, 0, 1, rel.tol = 1e-9)$value
  expect_equal(val, log(ref), tolerance = 1e-6)
})

test_that("k_cond: the default is 10 since 0.17.0 and k_cond = 6 is the fit of 0.16.0", {
  skip_on_cran()
  fx <- readRDS(test_path("fixtures", "nonsimplified_0.16.0.rds"))   # computed with 0.16.0
  d <- sim_ns(300)
  f6 <- tricopula(fm_ns, d, copula = "pcc", families = fam_ns, simplified = FALSE, k_cond = 6)
  expect_equal(coef(f6), fx$coef, tolerance = 1e-8)
  expect_equal(f6$logLik, fx$logLik, tolerance = 1e-10)
  # numerically differentiated covariance: platform differences of ~1e-5
  # relative (see test-regression.R); estimates and likelihood stay strict
  expect_equal(vcov(f6), fx$vcov, tolerance = 1e-3)
  expect_equal(unname(f6$sp), unname(fx$sp), tolerance = 1e-6)
  f10 <- tricopula(fm_ns, d, copula = "pcc", families = fam_ns, simplified = FALSE)
  expect_identical(f10$nonsimplified$k, 10L)
  expect_length(f10$index$ns, 9L)
  expect_equal(f10$convergence, 0L)
  # the same curve up to the larger basis: fitted dependence close to that of k = 6
  expect_lt(max(abs(f10$pair[, 3] - f6$pair[, 3])), 0.1)
  bt <- tc_bootstrap(f10, B = 2, seed = 1)
  expect_equal(colnames(bt$replicates)[1:21], names(coef(f10)))
})
