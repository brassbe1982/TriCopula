# Vine root (0.10.0)

sim_root <- function(n, fam, par, root, seed, df = 5) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n))
  U <- tc_rpcc(n, fam, par, df = df, root = root)
  d$y1 <- qnorm(U[, 1], 1 + 0.5 * d$x, 1.2)
  d$y2 <- qlnorm(U[, 2], 0.2 + 0.3 * d$x, 0.5)
  d$y3 <- qgamma(U[, 3], shape = 2.5, scale = exp(0.1 + 0.2 * d$x) / 2.5)
  d
}

test_that("vine density with each root matches VineCopula and its score is correct", {
  skip_on_cran()
  skip_if_not_installed("VineCopula")
  set.seed(1)
  U <- matrix(runif(240, 0.02, 0.98), 80, 3)
  W <- rep(list(matrix(1, 80, 1, dimnames = list(NULL, "(Intercept)"))), 3)
  fam <- c("frank", "student", "clayton")
  th <- c(3, atanh(0.4), log(1.2), log(6 - 2.05))
  for (r in 1:3) {
    cop <- structure("pcc", families = fam, root = r)
    ord <- TriCopula:::.tc_pcc_order(r)
    RVM <- VineCopula::C2RVine(ord, family = c(5, 2, 3), par = c(3, 0.4, 1.2),
                               par2 = c(0, 6, 0))
    ev <- TriCopula:::.tc_cop_eval(th, U, W, cop, deriv = TRUE)
    expect_equal(sum(ev$lc), VineCopula::RVineLogLik(U, RVM)$loglik, tolerance = 1e-9)
    g <- numDeriv::grad(function(x) sum(TriCopula:::.tc_cop_eval(x, U, W, cop)$lc), th)
    gu <- numDeriv::grad(function(x) {
      sum(TriCopula:::.tc_cop_eval(th, matrix(x, 80, 3), W, cop)$lc)
    }, as.vector(U))
    expect_equal(colSums(ev$score), g, tolerance = 1e-5)
    expect_equal(as.vector(ev$du), gu, tolerance = 1e-5)
  }
  expect_equal(TriCopula:::.tc_pcc_labels(1)$long, c("pcc12", "pcc13", "pcc23_given_1"))
  expect_equal(TriCopula:::.tc_pcc_labels(2)$short, c("12", "23", "13|2"))
  expect_equal(TriCopula:::.tc_pcc_labels(3)$long, c("pcc13", "pcc23", "pcc12_given_3"))
  # sampler: the root's pairs carry the dependence that was asked for
  set.seed(2)
  S <- tc_rpcc(2e4, c("clayton", "gumbel", "gaussian"), c(2, 2, 0), root = 3)
  expect_equal(cor(S[, 3], S[, 1], method = "kendall"), 0.5, tolerance = 0.02)
  expect_equal(cor(S[, 3], S[, 2], method = "kendall"), 0.5, tolerance = 0.02)
  expect_error(tc_rpcc(5, rep("gaussian", 3), c(0.1, 0.1, 0.1), root = 0), "root")
})

test_that("root = 1 is the default model; the root does not matter for a Gaussian vine", {
  skip_on_cran()
  d <- sim_mixed(400, seed = 151)
  f0 <- fit_mixed(d)
  g1 <- fit_mixed(d, "pcc", families = rep("gaussian", 3))
  g1b <- fit_mixed(d, "pcc", families = rep("gaussian", 3), root = 1)
  expect_identical(coef(g1b), coef(g1))
  expect_identical(g1b$vcov, g1$vcov)
  expect_equal(g1$root, 1L)
  for (r in 2:3) {
    gr <- fit_mixed(d, "pcc", families = rep("gaussian", 3), root = r)
    expect_equal(gr$logLik, f0$logLik, tolerance = 1e-6)
    expect_equal(unname(coef(gr))[1:9], unname(coef(f0))[1:9], tolerance = 1e-4)
  }
  expect_error(fit_mixed(d, "gaussian", root = 2), "only used")
  expect_error(fit_mixed(d, "pcc", families = rep("gaussian", 3), root = 4), "root must")
})

test_that("fit with root = 3: recovery, names, reordering equivalence, post-estimation", {
  skip_on_cran()
  fam <- c("clayton", "gumbel180", "frank")
  d <- sim_root(2000, fam, c(2, 1.8, 3), root = 3, seed = 152)
  fit <- fit_mixed(d, "pcc", families = fam, root = 3)
  expect_equal(fit$convergence, 0L)
  truth <- c(1, 0.5, log(1.2), 0.2, 0.3, log(0.5), 0.1, 0.2, log(2.5),
             log(2), log(0.8), 3)
  z <- (coef(fit) - truth) / sqrt(diag(vcov(fit)))
  expect_lt(max(abs(z)), 4)
  expect_lt(max(abs(fit$gradient)), 1e-2)
  expect_equal(names(coef(fit))[10:12],
               c("pcc13_clayton_log_theta", "pcc23_gumbel180_log_theta_minus_1",
                 "pcc12_given_3_frank_theta"))
  par <- truth * 1.02
  expect_equal(unname(tc_gradient(fit, par)),
               numDeriv::grad(function(p) tc_loglik(fit, p), par), tolerance = 1e-4)
  # the same model with the responses reordered so that y3 comes first
  alt <- tricopula(list(y3 ~ x, y1 ~ x, y2 ~ x), d,
                   margins = c("gamma", "gaussian", "lognormal"),
                   copula = "pcc", families = fam)
  expect_equal(fit$logLik, alt$logLik, tolerance = 1e-6)
  expect_equal(unname(coef(fit))[c(7:9, 1:6, 10:12)], unname(coef(alt)), tolerance = 1e-4)
  # a different root is a different model here
  f1 <- fit_mixed(d, "pcc", families = fam, hessian = FALSE)
  expect_gt(fit$logLik, f1$logLik + 5)
  # outputs carry the pair labels of the vine
  expect_named(tc_tau(fit), c("tau13", "tau23", "tau12|3"))
  expect_true(all(c("par13", "tau23", "par12|3") %in% rownames(tc_dependence(fit))))
  expect_named(predict(fit, data.frame(x = 0), type = "dependence"),
               c("par13", "par23", "par12|3", "tau13", "tau23", "tau12|3"))
  td <- tc_taildep(fit)
  expect_equal(td$pair, c("y3-y1", "y3-y2"))
  expect_gt(td$lower[1], 0)
  expect_output(print(fit), "c13, c23, c12\\|3")
  expect_lt(tc_diagnostics(fit)$max_abs_score_correlation_gap, 0.08)
  bt <- tc_bootstrap(fit, B = 4, seed = 4)
  expect_equal(bt$failures, 0L)
  # simulation from the fit reproduces the dependence of the data
  set.seed(153)
  s <- simulate(fit)
  expect_equal(cor(s$y1, s$y3, method = "kendall"), cor(d$y1, d$y3, method = "kendall"),
               tolerance = 0.05)
  expect_equal(cor(s$y1, s$y2, method = "kendall"), cor(d$y1, d$y2, method = "kendall"),
               tolerance = 0.05)
  skip_if_not_installed("mvtnorm")
  sn <- simulate(fit, newdata = data.frame(x = rep(0.3, 1e5)))
  th <- c(1, 1.2, 1)
  expect_equal(tc_jointprob(fit, th, data.frame(x = 0.3)),
               mean(sn$y1 <= th[1] & sn$y2 <= th[2] & sn$y3 <= th[3]), tolerance = 0.01)
  expect_equal(tc_jointprob(fit, th, data.frame(x = 0.3)),
               tc_jointprob(alt, th[c(3, 1, 2)], data.frame(x = 0.3)), tolerance = 1e-5)
})

test_that("tc_select with a given root and with root = \"auto\"", {
  skip_on_cran()
  fam <- c("clayton", "gumbel", "gaussian")
  d <- sim_root(1500, fam, c(3, 2.5, 0.05), root = 2, seed = 154)
  s2 <- tc_select(list(y1 ~ x, y2 ~ x, y3 ~ x), d, root = 2, hessian = FALSE)
  expect_equal(s2$root, 2L)
  expect_equal(s2$families[1:2], fam[1:2])
  expect_equal(unique(s2$table$pair), c("12", "23", "13|2"))
  expect_equal(s2$fit$root, 2L)
  sa <- tc_select(list(y1 ~ x, y2 ~ x, y3 ~ x), d, root = "auto", hessian = FALSE)
  expect_equal(nrow(sa$roots), 3L)
  expect_equal(sa$root, 2L)
  expect_equal(sa$roots$selected, c(FALSE, TRUE, FALSE))
  expect_equal(sa$roots$AIC[2], AIC(s2$fit), tolerance = 1e-6)
  expect_equal(sa$families, s2$families)
})
