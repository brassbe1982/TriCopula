vc_code <- c(gaussian = 1, student = 2, clayton = 3, gumbel = 4,
             clayton180 = 13, gumbel180 = 14, clayton90 = 23, gumbel90 = 24,
             clayton270 = 33, gumbel270 = 34)

test_that("pair-copula densities and h-functions match VineCopula", {
  skip_if_not_installed("VineCopula")
  set.seed(1)
  u <- runif(40, 0.02, 0.98); v <- runif(40, 0.02, 0.98)
  # the ten families of 0.7.0; Frank and Joe are tested in test-frank-joe.R
  for (fam in names(vc_code)) {
    base <- TriCopula:::.tc_pair_base(fam)
    eta <- 0.4
    nat <- TriCopula:::.tc_pair_natural(eta, fam)
    neg <- TriCopula:::.tc_pair_rot(fam) %in% c(90L, 270L)
    tau <- if (base == "student") log(5.5 - 2.05) else NULL
    p <- TriCopula:::.tc_pair(u, v, eta, fam, tau)
    par <- if (neg) -nat else nat
    par2 <- if (base == "student") 5.5 else 0
    expect_equal(p$l, log(VineCopula::BiCopPDF(u, v, vc_code[[fam]], par, par2)),
                 tolerance = 1e-8, info = fam)
    expect_equal(p$h, VineCopula::BiCopHfunc1(u, v, vc_code[[fam]], par, par2),
                 tolerance = 1e-8, info = fam)
    expect_equal(TriCopula:::.tc_pair_ktau(eta, fam),
                 VineCopula::BiCopPar2Tau(vc_code[[fam]], par, par2),
                 tolerance = 1e-8, info = fam)
    # inverse h-function
    expect_equal(TriCopula:::.tc_pair_hinv(p$h, u, eta, fam, tau), v,
                 tolerance = 1e-7, info = fam)
    # finite-difference derivatives against numDeriv
    d <- TriCopula:::.tc_pair_d(u, v, eta, fam, tau)
    f <- function(x, what) sum(TriCopula:::.tc_pair(x[1:40], x[41:80], x[81], fam, tau)[[what]])
    gl <- numDeriv::grad(f, c(u, v, eta), what = "l")
    gh <- numDeriv::grad(f, c(u, v, eta), what = "h")
    expect_equal(c(d$dl_u, d$dl_v, sum(d$dl_eta)), gl, tolerance = 1e-5, info = fam)
    expect_equal(c(d$dh_u, d$dh_v, sum(d$dh_eta)), gh, tolerance = 1e-5, info = fam)
  }
})

test_that("the C-vine reproduces the trivariate Gaussian and Student-t copulas", {
  set.seed(2)
  U <- matrix(runif(90, 0.02, 0.98), 30, 3)
  W <- rep(list(matrix(1, 30, 1, dimnames = list(NULL, "(Intercept)"))), 3)
  th <- TriCopula:::.tc_dep_inverse(R_true)
  g <- TriCopula:::.tc_pcc_eval(th, U, W, rep("gaussian", 3))$lc
  expect_equal(g, TriCopula:::.tc_logcop(U, R_true, Inf, "gaussian")$lc, tolerance = 1e-10)
  # trivariate t(nu) = vine with t pairs of df (nu, nu, nu + 1)
  nu <- 6
  tt <- TriCopula:::.tc_pcc_eval(c(th, log(nu - 2.05), log(nu - 2.05), log(nu + 1 - 2.05)),
                                 U, W, rep("student", 3))$lc
  expect_equal(tt, TriCopula:::.tc_logcop(U, R_true, nu, "student")$lc, tolerance = 1e-9)
})

test_that("the C-vine density integrates to one and the sampler matches it", {
  skip_on_cran()
  fam <- c("clayton", "gumbel180", "gumbel")
  W <- function(n) rep(list(matrix(1, n, 1, dimnames = list(NULL, "(Intercept)"))), 3)
  set.seed(3)
  U <- matrix(runif(6e5), ncol = 3)
  th <- c(log(1.5), log(0.8), log(0.4))
  expect_equal(mean(exp(TriCopula:::.tc_pcc_eval(th, U, W(nrow(U)), fam)$lc)), 1,
               tolerance = 0.05)
  S <- tc_rpcc(2e4, fam, c(1.5, 1.8, 1.4))
  expect_equal(cor(S[, 1], S[, 2], method = "kendall"), 1.5 / 3.5, tolerance = 0.02)
  expect_equal(cor(S[, 1], S[, 3], method = "kendall"), 1 - 1 / 1.8, tolerance = 0.02)
  # Rosenblatt transforms of simulated data are independent uniforms
  ro <- TriCopula:::.tc_pcc_eval(th, S, W(nrow(S)), fam)$rosenblatt
  expect_lt(max(abs(cor(qnorm(ro))[upper.tri(diag(3))])), 0.03)
  expect_gt(ks.test(ro[, 3], "punif")$p.value, 0.001)
  expect_error(tc_rpcc(5, c("clayton", "amh", "gumbel"), c(1, 1, 2)), "families")
  expect_error(tc_rpcc(5, c("clayton", "gaussian", "gumbel"), c(1, 0.2, 0.5)), "exceed 1")
})

sim_pcc <- function(n, fam, par, seed = 1, df = 5) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n))
  U <- tc_rpcc(n, fam, par, df)
  d$y1 <- qnorm(U[, 1], 1 + 0.5 * d$x, 1.2)
  d$y2 <- qlnorm(U[, 2], 0.2 + 0.3 * d$x, 0.5)
  d$y3 <- qgamma(U[, 3], shape = 2.5, scale = exp(0.1 + 0.2 * d$x) / 2.5)
  d
}

test_that("an all-Gaussian PCC fit equals the Gaussian-copula fit", {
  d <- sim_mixed(400, seed = 111)
  f0 <- fit_mixed(d)
  f1 <- fit_mixed(d, "pcc", families = rep("gaussian", 3))
  expect_equal(f1$logLik, f0$logLik, tolerance = 1e-6)
  expect_equal(unname(coef(f1)), unname(coef(f0)), tolerance = 1e-4)
  expect_equal(unname(sqrt(diag(vcov(f1)))), unname(sqrt(diag(vcov(f0)))), tolerance = 1e-3)
  expect_equal(names(coef(f1))[10:12],
               c("pcc12_gaussian_atanh_rho", "pcc13_gaussian_atanh_rho",
                 "pcc23_given_1_gaussian_atanh_rho"))
  expect_null(f1$correlation)
})

test_that("joint PCC fit recovers parameters and its gradient is correct", {
  skip_on_cran()
  fam <- c("clayton", "gumbel180", "student")
  d <- sim_pcc(2500, fam, c(2, 1.6, 0.4), seed = 112)
  fit <- fit_mixed(d, "pcc", families = fam)
  expect_equal(fit$convergence, 0L)
  expect_false(is.null(fit$vcov))
  truth <- c(1, 0.5, log(1.2), 0.2, 0.3, log(0.5), 0.1, 0.2, log(2.5),
             log(2), log(0.6), atanh(0.4), log(5 - 2.05))
  z <- (coef(fit) - truth) / sqrt(diag(vcov(fit)))
  expect_lt(max(abs(z)), 4)
  expect_gte(fit$logLik, fit$ifm$logLik - 1e-6)
  expect_lt(max(abs(fit$gradient)), 1e-2)
  par <- truth * 1.03
  expect_equal(unname(tc_gradient(fit, par)),
               numDeriv::grad(function(p) tc_loglik(fit, p), par), tolerance = 1e-4)
  expect_equal(names(coef(fit))[10:13],
               c("pcc12_clayton_log_theta", "pcc13_gumbel180_log_theta_minus_1",
                 "pcc23_given_1_student_atanh_rho", "pcc23_given_1_log_df_minus_2.05"))

  # the PCC beats the Gaussian copula on these data
  fg <- fit_mixed(d)
  cmp <- tc_compare(fg, fit)
  expect_lt(cmp$AIC[2], cmp$AIC[1])
  expect_match(cmp$copula[2], "pcc\\(clayton,gumbel180,student\\)")

  # post-estimation
  dep <- tc_dependence(fit)
  expect_true(all(c("par12", "tau13", "df23|1") %in% rownames(dep)))
  expect_equal(dep["tau12", "estimate"], dep["par12", "estimate"] / (dep["par12", "estimate"] + 2))
  expect_true(all(is.finite(dep$std.error)))
  expect_named(tc_tau(fit), c("tau12", "tau13", "tau23|1"))
  td <- tc_taildep(fit)
  expect_equal(td$lower[1], 2^(-1 / dep["par12", "estimate"]))
  expect_equal(td$upper[1], 0)
  expect_gt(td$lower[2], 0)       # survival Gumbel: lower-tail dependence
  pd <- predict(fit, data.frame(x = 0), type = "dependence")
  expect_equal(dim(pd), c(1L, 6L))
  dg <- tc_diagnostics(fit)
  expect_lt(dg$max_abs_score_correlation_gap, 0.08)
  expect_output(print(fit), "Pair copulas")
  expect_output(print(summary(fit)), "natural scale")

  skip_if_not_installed("mvtnorm")
  set.seed(113)
  nd <- data.frame(x = rep(0.3, 1e5))
  s <- simulate(fit, newdata = nd)
  th <- c(1, 1.2, 1)
  jp <- tc_jointprob(fit, th, data.frame(x = 0.3))
  expect_equal(jp, mean(s$y1 <= th[1] & s$y2 <= th[2] & s$y3 <= th[3]), tolerance = 0.01)
})

test_that("family selection recovers the generating families", {
  skip_on_cran()
  fam <- c("clayton", "gumbel", "gaussian")
  d <- sim_pcc(2000, fam, c(2, 1.8, 0.4), seed = 121)
  sel <- tc_select(list(y1 ~ x, y2 ~ x, y3 ~ x), d,
                   margins = c("gaussian", "lognormal", "gamma"))
  expect_equal(sel$families[1:2], fam[1:2])
  expect_true(sel$families[3] %in% c("gaussian", "student"))
  expect_s3_class(sel$fit, "tricopula")
  expect_equal(sel$fit$families, sel$families)
  expect_equal(nrow(sel$table), 18L)
  expect_equal(sum(sel$table$selected), 3L)
  expect_true("AIC" %in% names(sel$table))
  # a single pair fit agrees with VineCopula's ML estimate
  skip_if_not_installed("VineCopula")
  set.seed(122)
  U <- tc_rpcc(1500, c("clayton", "gaussian", "gaussian"), c(1.5, 0.3, 0.2))
  mine <- TriCopula:::.tc_fit_pair(U[, 1], U[, 2], "clayton")
  vc <- VineCopula::BiCopEst(U[, 1], U[, 2], family = 3, method = "mle")
  expect_equal(exp(mine$eta), vc$par, tolerance = 1e-3)
  expect_equal(mine$logLik, vc$logLik, tolerance = 1e-4)
  # negative dependence switches to the 90/270 rotations
  dn <- sim_pcc(600, c("clayton90", "gaussian", "gaussian"), c(2, 0.3, 0.2), seed = 123)
  sn <- tc_select(list(y1 ~ x, y2 ~ x, y3 ~ x), dn,
                  margins = c("gaussian", "lognormal", "gamma"),
                  criterion = "BIC", hessian = FALSE)
  expect_true(all(sn$table$family[sn$table$pair == "12"] %in%
                    c("gaussian", "student", "clayton90", "gumbel90", "clayton270", "gumbel270")))
  expect_equal(sn$families[1], "clayton90")
  expect_error(tc_select(list(y1 ~ x, y2 ~ x, y3 ~ x), dn, candidates = "amh"), "candidates")
})

test_that("PCC with dependence covariates, rotations, bootstrap and input checks", {
  skip_on_cran()
  set.seed(114)
  n <- 800
  d <- data.frame(x = rnorm(n))
  U <- TriCopula:::.tc_rpcc(n, list(log(1.5) + 0.4 * d$x, log(0.8), atanh(-0.3)),
                            c("clayton", "gumbel90", "gaussian"), list(NULL, NULL, NULL))
  d$y1 <- qnorm(U[, 1], 1 + 0.5 * d$x, 1.2)
  d$y2 <- qlnorm(U[, 2], 0.2 + 0.3 * d$x, 0.5)
  d$y3 <- qgamma(U[, 3], shape = 2.5, scale = exp(0.1 + 0.2 * d$x) / 2.5)
  fit <- fit_mixed(d, "pcc", families = c("clayton", "gumbel90", "gaussian"),
                   dependence = list(~ x, ~ 1, ~ 1))
  expect_true("pcc12_clayton_log_theta_x" %in% names(coef(fit)))
  expect_lt(tc_tau(fit)[["tau13"]], 0)
  b <- coef(fit); se <- sqrt(diag(vcov(fit)))
  expect_lt(abs(b["pcc12_clayton_log_theta_x"] - 0.4) / se["pcc12_clayton_log_theta_x"], 4)
  bt <- tc_bootstrap(fit, B = 6, seed = 1)
  expect_equal(bt$failures, 0L)
  ef <- tc_effects(fit, "x")
  expect_equal(ef$effect[1], unname(b["margin1_x"]), tolerance = 1e-6)
  expect_error(fit_mixed(d, "pcc"), "supply families")
  expect_error(fit_mixed(d, "pcc", families = c("clayton", "amh", "gumbel")), "supply families")
  expect_error(fit_mixed(d, "gaussian", families = rep("gaussian", 3)), "only used")
})
