# Frank and Joe pair copulas (0.10.0)

vc_fj <- c(frank = 5, joe = 6, joe180 = 16, joe90 = 26, joe270 = 36)

test_that("Frank and Joe densities, h-functions and tau match independent references", {
  skip_if_not_installed("VineCopula")
  set.seed(1)
  u <- runif(60, 0.02, 0.98); v <- runif(60, 0.02, 0.98)
  for (fam in names(vc_fj)) {
    etas <- if (fam == "frank") c(-34, -6, -0.5, 1e-3, 0.4, 3, 12, 34) else c(-3, -0.5, 0.4, 2.5)
    for (eta in etas) {
      nat <- TriCopula:::.tc_pair_natural(eta, fam)
      par <- if (TriCopula:::.tc_pair_rot(fam) %in% c(90L, 270L)) -nat else nat
      p <- TriCopula:::.tc_pair(u, v, eta, fam)
      info <- paste(fam, eta)
      expect_lt(max(abs(p$l - log(VineCopula::BiCopPDF(u, v, vc_fj[[fam]], par)))), 1e-9,
                label = info)
      expect_lt(max(abs(p$h - VineCopula::BiCopHfunc1(u, v, vc_fj[[fam]], par))), 1e-10,
                label = info)
      # VineCopula's Frank tau is itself approximate (see the next block)
      expect_lt(abs(TriCopula:::.tc_pair_ktau(eta, fam) -
                      VineCopula::BiCopPar2Tau(vc_fj[[fam]], par)),
                if (fam == "frank") 1e-3 else 1e-10, label = info)
    }
  }
  # Frank tau against the Debye integral computed by stats::integrate
  ref <- function(th) {
    a <- abs(th)
    I <- integrate(function(t) t / expm1(t), 0, a, rel.tol = 1e-13, abs.tol = 0)$value
    sign(th) * (1 - 4 / a * (1 - I / a))
  }
  for (th in c(-60, -30, -3, -0.5, 1e-3, 0.05, 0.0999, 0.1001, 0.4, 2.5, 34, 90)) {
    # absolute tolerance: the reference itself cancels near theta = 0
    expect_lt(abs(TriCopula:::.tc_frank_tau(th) - ref(th)), 1e-10, label = paste("frank tau", th))
  }
  expect_equal(TriCopula:::.tc_frank_tau(c(2, -2, 2)), c(ref(2), -ref(2), ref(2)),
               tolerance = 1e-10)
  # Joe tau at theta = 2 (removable singularity of the formula)
  expect_equal(TriCopula:::.tc_joe_tau(2), TriCopula:::.tc_joe_tau(2 + 1e-5), tolerance = 1e-5)
  # tail dependence
  for (fam in c("joe", "joe180", "joe90", "frank")) {
    eta <- 0.3
    nat <- TriCopula:::.tc_pair_natural(eta, fam)
    par <- if (TriCopula:::.tc_pair_rot(fam) %in% c(90L, 270L)) -nat else nat
    expect_equal(unname(TriCopula:::.tc_pair_taildep(eta, fam)),
                 unname(unlist(VineCopula::BiCopPar2TailDep(vc_fj[[fam]], par))),
                 tolerance = 1e-12, info = fam)
  }
})

test_that("Frank and Joe: inverse h, derivatives, tau inversion, normalization", {
  set.seed(2)
  u <- runif(40, 0.05, 0.95); v <- runif(40, 0.05, 0.95)
  for (fam in names(vc_fj)) {
    for (eta in if (fam == "frank") c(-4, 0.4, 5) else c(-1, 0.4)) {
      p <- TriCopula:::.tc_pair(u, v, eta, fam)
      expect_equal(TriCopula:::.tc_pair_hinv(p$h, u, eta, fam), v, tolerance = 1e-7,
                   info = paste(fam, eta))
      d <- TriCopula:::.tc_pair_d(u, v, eta, fam)
      f <- function(x, what) sum(TriCopula:::.tc_pair(x[1:40], x[41:80], x[81], fam)[[what]])
      gl <- numDeriv::grad(f, c(u, v, eta), what = "l")
      gh <- numDeriv::grad(f, c(u, v, eta), what = "h")
      expect_equal(c(d$dl_u, d$dl_v, sum(d$dl_eta)), gl, tolerance = 1e-5, info = fam)
      expect_equal(c(d$dh_u, d$dh_v, sum(d$dh_eta)), gh, tolerance = 1e-5, info = fam)
    }
    for (k in c(0.15, 0.5, 0.8)) {
      kt <- if (TriCopula:::.tc_pair_rot(fam) %in% c(90L, 270L)) -k else k
      expect_equal(TriCopula:::.tc_pair_ktau(TriCopula:::.tc_pair_eta_from_tau(kt, fam), fam),
                   kt, tolerance = 1e-6, info = fam)
    }
  }
  expect_equal(TriCopula:::.tc_pair_ktau(TriCopula:::.tc_pair_eta_from_tau(-0.4, "frank"), "frank"),
               -0.4, tolerance = 1e-6)
  # the density integrates to one, and h is a CDF in v
  g <- seq(0.0025, 0.9975, by = 0.005)
  G <- expand.grid(u = g, v = g)
  for (fam in c("frank", "joe", "joe90")) {
    eta <- if (fam == "frank") -3 else 0.2
    expect_equal(mean(exp(TriCopula:::.tc_pair(G$u, G$v, eta, fam)$l)), 1, tolerance = 5e-3,
                 info = fam)
  }
  expect_true(all(diff(TriCopula:::.tc_pair(0.3, g, 4, "frank")$h) > 0))
  expect_true(all(diff(TriCopula:::.tc_pair(0.3, g, 0.5, "joe")$h) > 0))
  # outside the numerical domain the density is rejected
  expect_true(all(is.infinite(TriCopula:::.tc_pair(u, v, 101, "frank")$l)))
  expect_true(all(is.finite(TriCopula:::.tc_pair(u, v, 99, "frank")$l)))
})

test_that("sampling and single-pair ML with Frank and Joe", {
  skip_on_cran()
  set.seed(3)
  S <- tc_rpcc(2e4, c("frank", "joe", "frank"), c(-5, 2, 3))
  expect_equal(cor(S[, 1], S[, 2], method = "kendall"), TriCopula:::.tc_frank_tau(-5),
               tolerance = 0.02)
  expect_equal(cor(S[, 1], S[, 3], method = "kendall"), TriCopula:::.tc_joe_tau(2),
               tolerance = 0.02)
  W <- rep(list(matrix(1, nrow(S), 1, dimnames = list(NULL, "(Intercept)"))), 3)
  ro <- TriCopula:::.tc_pcc_eval(c(-5, 0, 3), S, W, c("frank", "joe", "frank"))$rosenblatt
  expect_lt(max(abs(cor(qnorm(ro))[upper.tri(diag(3))])), 0.03)
  expect_error(tc_rpcc(5, c("frank", "joe", "gaussian"), c(0, 2, 0.1)), "nonzero")
  expect_error(tc_rpcc(5, c("frank", "joe", "gaussian"), c(1, 1, 0.1)), "Joe theta")
  expect_error(tc_rpcc(5, c("frank180", "joe", "gaussian"), c(1, 2, 0.1)), "families")
  skip_if_not_installed("VineCopula")
  U <- S[1:1500, ]
  mf <- TriCopula:::.tc_fit_pair(U[, 1], U[, 2], "frank")
  vf <- VineCopula::BiCopEst(U[, 1], U[, 2], family = 5, method = "mle")
  expect_equal(mf$eta, vf$par, tolerance = 1e-3)
  expect_equal(mf$logLik, vf$logLik, tolerance = 1e-4)
  mj <- TriCopula:::.tc_fit_pair(U[, 1], U[, 3], "joe")
  vj <- VineCopula::BiCopEst(U[, 1], U[, 3], family = 6, method = "mle")
  expect_equal(1 + exp(mj$eta), vj$par, tolerance = 1e-3)
  expect_equal(mj$logLik, vj$logLik, tolerance = 1e-4)
})

test_that("joint fit with Frank and Joe pairs", {
  skip_on_cran()
  fam <- c("frank", "joe180", "joe")
  set.seed(131)
  n <- 1500
  d <- data.frame(x = rnorm(n))
  U <- tc_rpcc(n, fam, c(-4, 2, 1.6))
  d$y1 <- qnorm(U[, 1], 1 + 0.5 * d$x, 1.2)
  d$y2 <- qlnorm(U[, 2], 0.2 + 0.3 * d$x, 0.5)
  d$y3 <- qgamma(U[, 3], shape = 2.5, scale = exp(0.1 + 0.2 * d$x) / 2.5)
  fit <- fit_mixed(d, "pcc", families = fam)
  expect_equal(fit$convergence, 0L)
  expect_false(is.null(fit$vcov))
  truth <- c(1, 0.5, log(1.2), 0.2, 0.3, log(0.5), 0.1, 0.2, log(2.5), -4, 0, log(0.6))
  z <- (coef(fit) - truth) / sqrt(diag(vcov(fit)))
  expect_lt(max(abs(z)), 4)
  expect_lt(max(abs(fit$gradient)), 1e-2)
  expect_gte(fit$logLik, fit$ifm$logLik - 1e-6)
  par <- truth + 0.03
  expect_equal(unname(tc_gradient(fit, par)),
               numDeriv::grad(function(p) tc_loglik(fit, p), par), tolerance = 1e-4)
  expect_equal(names(coef(fit))[10:12],
               c("pcc12_frank_theta", "pcc13_joe180_log_theta_minus_1",
                 "pcc23_given_1_joe_log_theta_minus_1"))
  dep <- tc_dependence(fit)
  expect_true(all(is.finite(dep$std.error)))
  expect_lt(dep["tau12", "estimate"], 0)
  expect_equal(dep["tau12", "estimate"],
               TriCopula:::.tc_frank_tau(dep["par12", "estimate"]), tolerance = 1e-8)
  td <- tc_taildep(fit)
  expect_equal(c(td$lower[1], td$upper[1]), c(0, 0))
  expect_gt(td$lower[2], 0); expect_equal(td$upper[2], 0)
  pd <- predict(fit, data.frame(x = 0), type = "dependence")
  expect_named(pd, c("par12", "par13", "par23|1", "tau12", "tau13", "tau23|1"))
  bt <- tc_bootstrap(fit, B = 4, seed = 2)
  expect_equal(bt$failures, 0L)
  skip_if_not_installed("mvtnorm")
  set.seed(132)
  s <- simulate(fit, newdata = data.frame(x = rep(0.3, 1e5)))
  th <- c(1, 1.2, 1)
  expect_equal(tc_jointprob(fit, th, data.frame(x = 0.3)),
               mean(s$y1 <= th[1] & s$y2 <= th[2] & s$y3 <= th[3]), tolerance = 0.01)
})

test_that("tc_select: default candidates are unchanged, \"all\" adds Frank and Joe", {
  skip_on_cran()
  set.seed(133)
  n <- 1200
  d <- data.frame(x = rnorm(n))
  U <- tc_rpcc(n, c("frank", "joe", "gaussian"), c(6, 2.2, 0.3))
  d$y1 <- qnorm(U[, 1], 1 + 0.5 * d$x, 1.2)
  d$y2 <- qlnorm(U[, 2], 0.2 + 0.3 * d$x, 0.5)
  d$y3 <- qgamma(U[, 3], shape = 2.5, scale = exp(0.1 + 0.2 * d$x) / 2.5)
  s0 <- tc_select(list(y1 ~ x, y2 ~ x, y3 ~ x), d, hessian = FALSE)
  expect_equal(nrow(s0$table), 18L)
  expect_false(any(s0$table$family %in% c("frank", "joe", "joe180")))
  expect_equal(s0$root, 1L)
  s1 <- tc_select(list(y1 ~ x, y2 ~ x, y3 ~ x), d, candidates = "all", hessian = FALSE)
  expect_equal(nrow(s1$table), 27L)
  expect_equal(s1$families[1:2], c("frank", "joe"))
  expect_gt(s1$fit$logLik, s0$fit$logLik)
  expect_error(tc_select(list(y1 ~ x, y2 ~ x, y3 ~ x), d, root = 4), "root")
})
