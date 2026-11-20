# Generalized Pareto, Dagum and Singh-Maddala margins (0.11.0)

me <- function(...) TriCopula:::.tc_margin_eval(...)

test_that("Dagum, Singh-Maddala and generalized Pareto match VGAM", {
  skip_if_not_installed("VGAM")
  set.seed(1)
  n <- 50
  X <- cbind("(Intercept)" = 1, x = rnorm(n))
  beta <- c(0.3, 0.4)
  eta <- drop(X %*% beta)
  y <- rexp(n) * 2 + 0.05
  extra <- log(1.7) + 0.1 * X[, 2]; shape <- log(0.6) - 0.2 * X[, 2]
  ev <- me(y, X, "dagum", beta, extra, shape = shape)
  expect_equal(ev$logf, VGAM::ddagum(y, exp(eta), exp(extra), exp(shape), log = TRUE),
               tolerance = 1e-12)
  expect_equal(ev$u, VGAM::pdagum(y, exp(eta), exp(extra), exp(shape)), tolerance = 1e-12)
  expect_equal(TriCopula:::.tc_margin_quantile(ev$u, eta, extra, "dagum", shape), y,
               tolerance = 1e-10)
  ev <- me(y, X, "singhmaddala", beta, extra, shape = shape)
  expect_equal(ev$logf, VGAM::dsinmad(y, exp(eta), exp(extra), exp(shape), log = TRUE),
               tolerance = 1e-12)
  expect_equal(ev$u, VGAM::psinmad(y, exp(eta), exp(extra), exp(shape)), tolerance = 1e-12)
  expect_equal(TriCopula:::.tc_margin_quantile(ev$u, eta, extra, "singhmaddala", shape), y,
               tolerance = 1e-10)
  for (xi in c(-0.4, -0.05, 0.3, 1.5)) {
    ex <- log(xi + 0.5)
    yy <- if (xi < 0) pmin(y, 0.95 * exp(eta) / abs(xi)) else y
    ev <- me(yy, X, "gpareto", beta, ex)
    expect_equal(ev$logf, VGAM::dgpd(yy, 0, exp(eta), xi, log = TRUE), tolerance = 1e-12,
                 info = xi)
    expect_equal(ev$u, VGAM::pgpd(yy, 0, exp(eta), xi), tolerance = 1e-12, info = xi)
    expect_equal(TriCopula:::.tc_margin_quantile(ev$u, eta, ex, "gpareto"), yy,
                 tolerance = 1e-10, info = xi)
  }
  # xi = 0 is the exponential distribution; vector and scalar predictors agree
  ev0 <- me(y, X, "gpareto", beta, log(0.5))
  expect_equal(ev0$logf, dexp(y, 1 / exp(eta), log = TRUE), tolerance = 1e-12)
  expect_equal(ev0$u, pexp(y, 1 / exp(eta)), tolerance = 1e-12)
  expect_equal(me(y, X, "gpareto", beta, rep(log(0.8), n))$logf,
               me(y, X, "gpareto", beta, log(0.8))$logf)
  # outside the support (xi < 0) the density is zero and the CDF one
  expect_identical(me(c(1, 30), X[1:2, ], "gpareto", c(0, 0), log(0.1))$logf[2], -Inf)
  expect_equal(TriCopula:::.tc_margin_cdf(30, 0, log(0.1), "gpareto"), 1)
})

test_that("new margins: normalization, means, scores and special cases", {
  for (m in c("dagum", "singhmaddala")) {
    f <- function(t) exp(me(t, matrix(1, length(t), 1), m, log(1.3), log(2.5),
                            shape = log(0.8))$logf)
    expect_equal(integrate(f, 0, Inf)$value, 1, tolerance = 1e-6)
    expect_equal(integrate(function(t) t * f(t), 0, Inf, rel.tol = 1e-10)$value,
                 TriCopula:::.tc_margin_mean(log(1.3), log(2.5), m, log(0.8)),
                 tolerance = 1e-6)
  }
  g <- function(t) exp(me(t, matrix(1, length(t), 1), "gpareto", log(1.3), log(0.7))$logf)
  expect_equal(integrate(g, 0, Inf)$value, 1, tolerance = 1e-6)
  expect_equal(integrate(function(t) t * g(t), 0, Inf, rel.tol = 1e-10)$value,
               TriCopula:::.tc_margin_mean(log(1.3), log(0.7), "gpareto"), tolerance = 1e-5)
  expect_equal(TriCopula:::.tc_margin_mean(0, log(0.9), "dagum", 0), Inf)        # a <= 1
  expect_equal(TriCopula:::.tc_margin_mean(0, log(1.6), "gpareto"), Inf)         # xi >= 1
  # analytic scores against numerical differentiation
  set.seed(2)
  n <- 30
  X <- cbind("(Intercept)" = 1, x = rnorm(n))
  y <- rexp(n) + 0.1
  chk <- function(m, beta, extra, shape = NULL) {
    ev <- me(y, X, m, beta, extra, deriv = TRUE, shape = shape)
    for (q in c("logf", "u")) {
      pre <- if (q == "logf") "dl_" else "du_"
      f <- function(x) sum(me(y, X, m, beta + c(x[1], 0), extra + x[2],
                              shape = if (!is.null(shape)) shape + x[3])[[q]])
      gr <- numDeriv::grad(f, c(0, 0, 0))
      expect_equal(sum(ev[[paste0(pre, "eta")]]), gr[1], tolerance = 1e-6, info = paste(m, q))
      expect_equal(sum(ev[[paste0(pre, "extra")]]), gr[2], tolerance = 1e-6, info = paste(m, q))
      if (!is.null(shape)) {
        expect_equal(sum(ev[[paste0(pre, "shape")]]), gr[3], tolerance = 1e-6, info = paste(m, q))
      }
    }
  }
  chk("dagum", c(0.2, 0.3), log(2) + 0.1 * X[, 2], log(0.7) - 0.2 * X[, 2])
  chk("singhmaddala", c(0.2, 0.3), log(2) + 0.1 * X[, 2], log(1.4) - 0.2 * X[, 2])
  chk("gpareto", c(0.2, 0.3), log(0.8))
  chk("gpareto", c(1.5, 0.1), log(0.3))            # xi = -0.2, all y inside the support
  chk("gpareto", c(0.2, 0.3), log(0.5))            # xi = 0 (series branch)
  # Dagum with p = 1 and Singh-Maddala with q = 1 are both the log-logistic
  d1 <- me(y, X, "dagum", c(0.2, 0.3), log(2), shape = 0)
  s1 <- me(y, X, "singhmaddala", c(0.2, 0.3), log(2), shape = 0)
  expect_equal(d1$logf, s1$logf, tolerance = 1e-12)
  expect_equal(d1$u, s1$u, tolerance = 1e-12)
  expect_equal(d1$u, plogis(2 * (log(y) - drop(X %*% c(0.2, 0.3)))), tolerance = 1e-12)
  # Singh-Maddala with a = 1 is the Lomax: generalized Pareto with xi = 1 / q
  q <- 2.5
  sm <- me(y, X, "singhmaddala", c(0.2, 0.3), 0, shape = log(q))
  gp <- me(y, X, "gpareto", c(0.2 - log(q), 0.3), log(1 / q + 0.5))
  expect_equal(sm$logf, gp$logf, tolerance = 1e-12)
  expect_equal(sm$u, gp$u, tolerance = 1e-12)
  # extreme parameters (found in run 11): the Dagum log-density must not
  # degrade to rounding noise when p is huge and the scale tiny
  ex <- me(y, X, "dagum", c(-35, -3), 4.2, deriv = TRUE, shape = 64)
  expect_true(all(is.finite(ex$logf)))
  expect_lt(sum(ex$logf), -1e4)
  expect_true(all(is.finite(ex$dl_eta)))
  big <- me(y, X, "dagum", c(0.2, 0.3), log(2), shape = log(400))
  expect_equal(big$logf, log(2) + log(400) - log(y) +
                 400 * 2 * (log(y) - drop(X %*% c(0.2, 0.3))) -
                 401 * log1p(exp(2 * (log(y) - drop(X %*% c(0.2, 0.3))))), tolerance = 1e-9)
  expect_error(me(y, X, "dagum", c(0, 0), 0), "shape predictor")
  expect_error(TriCopula:::.tc_check_support(c(1, -1), "gpareto"), "strictly positive")
})

sim_m3 <- function(n, seed, copula = "gaussian") {
  set.seed(seed)
  d <- data.frame(x = rnorm(n))
  U <- tc_rcopula(n, R_true, copula, 5)
  q <- TriCopula:::.tc_margin_quantile
  d$y1 <- q(U[, 1], 0.5 + 0.3 * d$x, log(3), "dagum", log(0.7))
  d$y2 <- q(U[, 2], 0.2 + 0.3 * d$x, log(2.5), "singhmaddala", log(1.5))
  d$y3 <- q(U[, 3], 0.1 + 0.2 * d$x, log(0.2 + 0.5), "gpareto")
  d
}
m3 <- c("dagum", "singhmaddala", "gpareto")
truth_m3 <- c(0.5, 0.3, log(3), log(0.7), 0.2, 0.3, log(2.5), log(1.5),
              0.1, 0.2, log(0.7))

test_that("shape = ~ 1 changes nothing for two-parameter margins", {
  d <- sim_mixed(250, seed = 161)
  f0 <- fit_mixed(d, hessian = FALSE)
  f1 <- fit_mixed(d, hessian = FALSE, shape = ~ 1)
  expect_identical(coef(f1), coef(f0))
  expect_identical(f1$logLik, f0$logLik)
  expect_null(f0$shape_design)
  expect_equal(lengths(f0$index$shapel), c(0L, 0L, 0L))
  expect_identical(f0$index$margin, list(1:3, 4:6, 7:9))
  expect_error(fit_mixed(d, shape = ~ x), "three-parameter")
})

test_that("joint fit with three-parameter margins: recovery and post-estimation", {
  skip_on_cran()
  d <- sim_m3(3000, 162)
  fit <- tricopula(list(y1 ~ x, y2 ~ x, y3 ~ x), d, margins = m3)
  expect_equal(fit$convergence, 0L)
  expect_false(is.null(fit$vcov))
  truth <- c(truth_m3, TriCopula:::.tc_dep_inverse(R_true))
  z <- (coef(fit) - truth) / sqrt(diag(vcov(fit)))
  expect_lt(max(abs(z)), 4)
  expect_lt(max(abs(fit$gradient)), 1e-2)
  expect_gte(fit$logLik, fit$ifm$logLik - 1e-6)
  expect_equal(names(coef(fit))[1:11],
               c("margin1_(Intercept)", "margin1_x", "margin1_log_a", "margin1_log_p",
                 "margin2_(Intercept)", "margin2_x", "margin2_log_a", "margin2_log_q",
                 "margin3_(Intercept)", "margin3_x", "margin3_log_shape_plus_0.5"))
  expect_equal(fit$index$margin, list(1:4, 5:8, 9:11))
  par <- truth * 1.03
  expect_equal(unname(tc_gradient(fit, par)),
               numDeriv::grad(function(p) tc_loglik(fit, p), par), tolerance = 1e-5)
  # log-likelihood against VGAM densities and the Gaussian copula from mvtnorm
  skip_if_not_installed("VGAM"); skip_if_not_installed("mvtnorm")
  b <- cbind(exp(0.5 + 0.3 * d$x), exp(0.2 + 0.3 * d$x), exp(0.1 + 0.2 * d$x))
  lm <- VGAM::ddagum(d$y1, b[, 1], 3, 0.7, log = TRUE) +
    VGAM::dsinmad(d$y2, b[, 2], 2.5, 1.5, log = TRUE) +
    VGAM::dgpd(d$y3, 0, b[, 3], 0.2, log = TRUE)
  zz <- qnorm(cbind(VGAM::pdagum(d$y1, b[, 1], 3, 0.7), VGAM::psinmad(d$y2, b[, 2], 2.5, 1.5),
                    VGAM::pgpd(d$y3, 0, b[, 3], 0.2)))
  lc <- mvtnorm::dmvnorm(zz, sigma = R_true, log = TRUE) - rowSums(dnorm(zz, log = TRUE))
  expect_equal(tc_loglik(fit, truth), sum(lm + lc), tolerance = 1e-9)
  # post-estimation
  dep <- tc_dependence(fit)
  expect_true(all(c("margin1_a", "margin1_p", "margin2_q", "margin3_shape_plus_0.5") %in%
                    rownames(dep)))
  expect_equal(dep["margin1_p", "estimate"], exp(coef(fit)[["margin1_log_p"]]))
  expect_true(all(is.finite(dep$std.error)))
  pp <- predict(fit, data.frame(x = 0))
  expect_named(pp$margin1, c("scale", "a", "p"))
  expect_named(pp$margin2, c("scale", "a", "q"))
  expect_equal(pp$margin3$shape, exp(coef(fit)[[11]]) - 0.5)
  pm <- predict(fit, data.frame(x = 0), type = "mean")
  set.seed(163)
  s <- simulate(fit, newdata = data.frame(x = rep(0, 2e5)))
  expect_equal(unname(colMeans(s)), unname(unlist(pm)), tolerance = 0.03)
  th <- c(1.5, 1, 1)
  expect_equal(tc_jointprob(fit, th, data.frame(x = 0)),
               mean(s$y1 <= th[1] & s$y2 <= th[2] & s$y3 <= th[3]), tolerance = 0.01)
  ef <- tc_effects(fit, "x")
  expect_true(all(is.finite(ef$std.error)))
  bt <- tc_bootstrap(fit, B = 4, seed = 5)
  expect_equal(bt$failures, 0L)
  expect_s3_class(tc_diagnostics(fit), "tc_diagnostics")
  expect_output(print(summary(fit)), "margin1_log_p")
})

test_that("shape predictor with covariates; three-parameter margins in a PCC", {
  skip_on_cran()
  set.seed(164)
  n <- 2500
  d <- data.frame(x = rnorm(n), g = rbinom(n, 1, 0.5))
  U <- tc_rpcc(n, c("clayton", "frank", "gaussian"), c(1.5, 3, 0.3))
  q <- TriCopula:::.tc_margin_quantile
  d$y1 <- q(U[, 1], 0.5 + 0.3 * d$x, log(3), "dagum", log(0.6) + 0.7 * d$g)
  d$y2 <- q(U[, 2], 0.2 + 0.3 * d$x, log(2.5), "singhmaddala", log(1.5))
  d$y3 <- qlnorm(U[, 3], 0.1 + 0.2 * d$x, 0.5)
  mg <- c("dagum", "singhmaddala", "lognormal")
  fam <- c("clayton", "frank", "gaussian")
  fit <- tricopula(list(y1 ~ x, y2 ~ x, y3 ~ x), d, margins = mg, copula = "pcc",
                   families = fam, shape = list(~ g, ~ 1, ~ 1))
  expect_equal(fit$convergence, 0L)
  expect_equal(names(coef(fit))[4:5], c("margin1_log_p_(Intercept)", "margin1_log_p_g"))
  truth <- c(0.5, 0.3, log(3), log(0.6), 0.7, 0.2, 0.3, log(2.5), log(1.5),
             0.1, 0.2, log(0.5), log(1.5), 3, atanh(0.3))
  z <- (coef(fit) - truth) / sqrt(diag(vcov(fit)))
  expect_lt(max(abs(z)), 4)
  expect_equal(unname(tc_gradient(fit, truth)),
               numDeriv::grad(function(p) tc_loglik(fit, p), truth), tolerance = 1e-4)
  # a zero slope gives the constant-shape likelihood
  f0 <- tricopula(list(y1 ~ x, y2 ~ x, y3 ~ x), d, margins = mg, copula = "pcc",
                  families = fam, hessian = FALSE)
  p0 <- coef(f0)
  expect_equal(tc_loglik(fit, c(p0[1:4], 0, p0[5:14])), f0$logLik, tolerance = 1e-10)
  expect_lt(tc_lrtest(f0, fit)$p.value, 0.01)
  pp <- predict(fit, data.frame(x = 0, g = 0:1))
  expect_equal(pp$margin1$p, exp(coef(fit)[[4]] + coef(fit)[[5]] * 0:1), tolerance = 1e-12)
  expect_equal(dim(simulate(fit, seed = 1)), c(n, 3L))
  sel <- tc_select(list(y1 ~ x, y2 ~ x, y3 ~ x), d, margins = mg,
                   shape = list(~ g, ~ 1, ~ 1), hessian = FALSE)
  expect_equal(sel$families[1], "clayton")
  expect_true("margin1_log_p_g" %in% names(coef(sel$fit)))
})
