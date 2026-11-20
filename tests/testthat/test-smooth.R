# Regression splines and penalized smooth terms (0.12.0)

sim_smooth <- function(n, seed, f1 = function(x) sin(2 * x)) {
  set.seed(seed)
  d <- data.frame(x = runif(n, -2, 2), z = rnorm(n))
  U <- tc_rcopula(n, R_true)
  d$y1 <- qnorm(U[, 1], 1 + f1(d$x) + 0.3 * d$z, 0.7)
  d$y2 <- qlnorm(U[, 2], 0.2 + 0.3 * d$x, 0.5)
  d$y3 <- qgamma(U[, 3], shape = 2.5, scale = exp(0.1 + 0.2 * d$x) / 2.5)
  d
}

test_that("stage A: unpenalized regression splines work in every formula", {
  d <- sim_smooth(600, 171)
  fit <- tricopula(list(y1 ~ z + splines::ns(x, 4), y2 ~ x, y3 ~ x), d,
                   scale = list(~ splines::ns(x, 2), ~ 1, ~ 1),
                   dependence = list(~ splines::bs(x, 4), ~ 1, ~ 1))
  expect_equal(fit$convergence, 0L)
  expect_null(fit$penalties); expect_null(fit$edf); expect_null(fit$sp)
  expect_equal(attr(logLik(fit), "df"), length(coef(fit)))
  # prediction on new data reuses the knots of the fit
  nd <- d[c(5, 50, 300), ]
  expect_equal(predict(fit, nd, type = "mean"), predict(fit, type = "mean")[c(5, 50, 300), ],
               ignore_attr = TRUE, tolerance = 1e-12)
  expect_equal(predict(fit, nd, type = "dependence")$rho12, unname(fit$rho[c(5, 50, 300), 1]),
               tolerance = 1e-12)
  expect_equal(predict(fit, nd)$margin1$scale,
               unname(fit$fitted_margins[[1]]$pars$sigma[c(5, 50, 300)]), tolerance = 1e-12)
  # the spline recovers the curve
  g <- data.frame(x = seq(-1.8, 1.8, length.out = 30), z = 0)
  expect_lt(max(abs(predict(fit, g, type = "mean")$y1 - (1 + sin(2 * g$x)))), 0.35)
  ef <- tc_effects(fit, "x", newdata = d[1:50, ])
  expect_true(all(is.finite(ef$effect)))
  expect_equal(dim(simulate(fit, seed = 1, newdata = nd)), c(3L, 3L))
})

test_that("formulas without s() keep the unpenalized path; input checks", {
  d <- sim_smooth(300, 172)
  f0 <- tricopula(list(y1 ~ z + x, y2 ~ x, y3 ~ x), d, hessian = FALSE)
  expect_null(f0$penalties); expect_null(f0$sp); expect_null(f0$edf)
  expect_identical(attr(logLik(f0), "df"), f0$parameter_count)
  expect_null(summary(f0)$smooths)
  expect_error(tricopula(list(y1 ~ z + x, y2 ~ x, y3 ~ x), d, sp = 1), "only used with s")
  expect_error(tc_smooth(f0), "no s\\(\\) terms")
  skip_if_not_installed("mgcv")
  expect_error(tricopula(list(y1 ~ t2(x, z), y2 ~ x, y3 ~ x), d), "t2")
  expect_error(tricopula(list(y1 ~ s(x, by = nope), y2 ~ x, y3 ~ x), d), "columns of the data")
  expect_error(vcov(f0, type = "corrected"), "no penalized")
  expect_error(tricopula(list(y1 ~ s(x), y2 ~ x, y3 ~ x), d, sp = c(1, 2)), "sp must be")
})

test_that("penalized fit: stationarity, fixed point, limits and agreement with mgcv", {
  skip_on_cran()
  skip_if_not_installed("mgcv")
  d <- sim_smooth(800, 173)
  fit <- tricopula(list(y1 ~ z + s(x), y2 ~ x, y3 ~ x), d)
  expect_equal(fit$convergence, 0L)
  expect_true(fit$sp_converged)
  expect_equal(names(fit$sp), "margin1:s(x)")
  expect_equal(names(coef(fit))[3:4], c("margin1_s(x).1", "margin1_s(x).2"))
  ii <- fit$penalties[[1]]$ind
  expect_equal(ii, 3:11)
  # independent rebuild of basis and penalty
  sm <- mgcv::smoothCon(mgcv::s(x), d, absorb.cons = TRUE, scale.penalty = TRUE)[[1]]
  expect_equal(unname(fit$design[[1]][, ii]), unname(sm$X), tolerance = 1e-12)
  S <- sm$S[[1]]
  b <- coef(fit)[ii]
  # stationarity of the penalized log-likelihood
  g <- tc_gradient(fit)
  expect_lt(max(abs(g[ii] - fit$sp * drop(S %*% b))), 1e-4)
  expect_lt(max(abs(g[-ii])), 1e-4)
  expect_lt(max(abs(fit$gradient)), 1e-4)
  expect_equal(unname(fit$penalized_logLik), fit$logLik - 0.5 * fit$sp[[1]] * sum(b * (S %*% b)),
               tolerance = 1e-10)
  # Bayesian covariance, edf, and the Fellner-Schall fixed point
  Hu <- -numDeriv::jacobian(function(p) tc_gradient(fit, p), coef(fit))
  Sl <- matrix(0, length(g), length(g)); Sl[ii, ii] <- fit$sp * S
  Vp <- solve((Hu + t(Hu)) / 2 + Sl)
  expect_equal(unname(vcov(fit)), Vp, tolerance = 1e-5)
  expect_equal(sum(fit$edf), sum(diag(Vp %*% Hu)), tolerance = 1e-4)
  upd <- (sm$rank / fit$sp - sum(diag(Vp[ii, ii] %*% S))) / sum(b * (S %*% b))
  expect_equal(unname(upd), 1, tolerance = 0.02)
  expect_equal(attr(logLik(fit), "df"), sum(fit$edf))
  expect_lt(sum(fit$edf), length(coef(fit)))
  # agreement with mgcv (REML) for the Gaussian margin
  gm <- mgcv::gam(y1 ~ z + s(x), data = d, method = "REML")
  cur <- tc_smooth(fit, 1, n = 60)
  pg <- predict(gm, data.frame(x = cur$x, z = 0), type = "terms", se.fit = TRUE)
  expect_lt(max(abs(cur$fit - pg$fit[, "s(x)"])), 0.15)
  expect_equal(sum(fit$edf[ii]), sum(gm$edf[-(1:2)]), tolerance = 0.15)
  expect_equal(fit$sp[[1]] * exp(2 * coef(fit)[["margin1_log_scale"]]), gm$sp[[1]],
               tolerance = 0.3)
  expect_true(all(cur$se / pg$se.fit[, "s(x)"] > 0.7 & cur$se / pg$se.fit[, "s(x)"] < 1.3))
  expect_gt(mean(cur$lower <= sin(2 * cur$x) - mean(sin(2 * d$x)) &
                   sin(2 * cur$x) - mean(sin(2 * d$x)) <= cur$upper), 0.8)
  # a very large smoothing parameter gives the linear model (null space of
  # the thin-plate penalty), sp = 0 the unpenalized spline
  lin <- tricopula(list(y1 ~ z + x, y2 ~ x, y3 ~ x), d, hessian = FALSE)
  big <- tricopula(list(y1 ~ z + s(x), y2 ~ x, y3 ~ x), d, sp = 1e7, hessian = FALSE)
  expect_equal(big$logLik, lin$logLik, tolerance = 1e-5)
  expect_equal(sum(big$edf), length(coef(lin)), tolerance = 1e-3)
  zero <- tricopula(list(y1 ~ z + s(x), y2 ~ x, y3 ~ x), d, sp = 0, hessian = FALSE)
  expect_equal(sum(zero$edf), length(coef(zero)), tolerance = 1e-6)
  expect_gt(zero$logLik, fit$logLik)
  expect_gt(fit$logLik, big$logLik)
  expect_lt(AIC(fit), AIC(big))
  # refit with the selected sp held fixed reproduces the fit
  re <- tricopula(list(y1 ~ z + s(x), y2 ~ x, y3 ~ x), d, sp = fit$sp, hessian = FALSE)
  expect_equal(coef(re), coef(fit), tolerance = 1e-4)
  # summary and prediction
  s <- summary(fit)
  expect_equal(s$smooths$term, "margin1:s(x)")
  expect_lt(s$smooths$p.value, 1e-6)
  expect_output(print(s), "Penalized smooth terms")
  nd <- d[c(3, 30, 300), ]
  expect_equal(predict(fit, nd, type = "mean"), predict(fit, type = "mean")[c(3, 30, 300), ],
               ignore_attr = TRUE, tolerance = 1e-10)
  expect_error(tc_smooth(fit, "nope"), "which must be")
  bt <- tc_bootstrap(fit, B = 3, seed = 6)
  expect_equal(bt$failures, 0L)
  # nonparametric bootstrap keeps the spline basis of the fit
  bn <- tc_bootstrap(fit, B = 3, type = "nonparametric", seed = 7)
  expect_equal(bn$failures, 0L)
  expect_equal(rownames(bn$coefficients), names(coef(fit)))
  # reusing the basis on the same data reproduces the design and the fit
  rb <- tricopula(list(y1 ~ z + s(x), y2 ~ x, y3 ~ x), d, sp = fit$sp, hessian = FALSE,
                  smooth_basis = fit$smooths)
  expect_equal(rb$design[[1]], fit$design[[1]], tolerance = 1e-8)
  expect_equal(coef(rb), coef(fit), tolerance = 1e-4)
  # covariance corrected for smoothing-parameter uncertainty
  Vc <- vcov(fit, type = "corrected")
  expect_equal(dim(Vc), dim(vcov(fit)))
  dV <- Vc - vcov(fit)
  expect_gt(min(eigen((dV + t(dV)) / 2, symmetric = TRUE, only.values = TRUE)$values), -1e-8)
  expect_true(all(diag(Vc)[ii] >= diag(vcov(fit))[ii] - 1e-12))
  # dtheta / drho used by the correction, against a refit at a nearby sp
  hh <- 0.05
  f2 <- tricopula(list(y1 ~ z + s(x), y2 ~ x, y3 ~ x), d, sp = fit$sp * exp(hh),
                  hessian = FALSE, start = coef(fit))
  J <- -fit$sp[[1]] * drop(vcov(fit)[, ii] %*% (S %*% b))
  expect_equal(unname((coef(f2) - coef(fit)) / hh)[ii], unname(J[ii]), tolerance = 0.1)
  cc <- tc_smooth(fit, 1, n = 20, vcov = Vc)
  cb <- tc_smooth(fit, 1, n = 20)
  expect_true(all(cc$se >= cb$se - 1e-12))
  expect_lt(max(cc$se / cb$se), 1.5)
})

test_that("smooths in the scale, dependence and df predictors", {
  skip_on_cran()
  skip_if_not_installed("mgcv")
  set.seed(174)
  n <- 2000
  d <- data.frame(x = runif(n, -2, 2))
  r12 <- tanh(0.6 + 0.6 * cos(1.5 * d$x))
  sig <- exp(-0.3 + 0.4 * sin(1.5 * d$x))
  U <- TriCopula:::.tc_rcopula_vec(n, r12, rep(0.3, n), rep(0.3, n), "gaussian", Inf)
  d$y1 <- qnorm(U[, 1], 1 + 0.5 * d$x, sig)
  d$y2 <- qlnorm(U[, 2], 0.2 + 0.3 * d$x, 0.5)
  d$y3 <- qgamma(U[, 3], shape = 2.5, scale = exp(0.1 + 0.2 * d$x) / 2.5)
  fit <- tricopula(list(y1 ~ x, y2 ~ x, y3 ~ x), d,
                   scale = list(~ s(x, k = 8), ~ 1, ~ 1),
                   dependence = list(~ s(x, k = 8), ~ 1, ~ 1))
  expect_equal(fit$convergence, 0L)
  expect_true(fit$sp_converged)
  expect_equal(names(fit$sp), c("margin1_log_scale:s(x)", "dep12:s(x)"))
  expect_true("atanh_partial_rho12_s(x).1" %in% names(coef(fit)))
  expect_true("margin1_log_scale_s(x).1" %in% names(coef(fit)))
  g <- data.frame(x = seq(-1.8, 1.8, length.out = 25))
  pd <- predict(fit, g, type = "dependence")
  expect_lt(max(abs(pd$rho12 - tanh(0.6 + 0.6 * cos(1.5 * g$x)))), 0.12)
  expect_lt(max(abs(log(predict(fit, g)$margin1$scale) - (-0.3 + 0.4 * sin(1.5 * g$x)))), 0.15)
  expect_equal(predict(fit, d[1:4, , drop = FALSE], type = "dependence")$rho12,
               unname(fit$rho[1:4, 1]), tolerance = 1e-10)
  tab <- summary(fit)$smooths
  expect_true(all(tab$edf > 1.5 & tab$edf < 7))
  expect_true(all(tab$p.value < 0.01))
  # the smooth model beats constant scale and dependence by AIC
  f0 <- tricopula(list(y1 ~ x, y2 ~ x, y3 ~ x), d, hessian = FALSE)
  expect_lt(AIC(fit), AIC(f0))
  cur <- tc_smooth(fit, "dep12:s(x)", n = 20)
  expect_equal(nrow(cur), 20L)
  expect_true(all(cur$se > 0))
  expect_equal(dim(simulate(fit, seed = 1)), c(n, 3L))
  skip_if_not_installed("mvtnorm")
  jp <- tc_jointprob(fit, c(1, 1.2, 1), data.frame(x = c(-1, 1)))
  expect_true(all(jp > 0 & jp < 1))

  # a flat truth is shrunk towards the constant: edf of the smooth near 1
  d2 <- sim_mixed(1200, "student", seed = 175)
  ft <- fit_mixed(d2, "student", df = ~ s(x, k = 6))
  expect_equal(names(ft$sp), "df:s(x)")
  expect_lt(summary(ft)$smooths$edf, 2.5)
  expect_equal(predict(ft, d2[1:3, ], type = "dependence")$df, ft$df_i[1:3], tolerance = 1e-10)
})

test_that("tensor products, by variables and smooths with several penalties", {
  skip_on_cran()
  skip_if_not_installed("mgcv")
  set.seed(176)
  n <- 900
  d <- data.frame(x = runif(n, -2, 2), z = runif(n, -1, 1),
                  g = factor(sample(c("a", "b"), n, TRUE)), w = rbinom(n, 1, 0.5))
  U <- tc_rcopula(n, R_true)
  f2 <- function(x, z) sin(1.5 * x) * z
  d$y1 <- qnorm(U[, 1], 1 + f2(d$x, d$z), 0.5)
  d$y2 <- qlnorm(U[, 2], 0.2 + ifelse(d$g == "a", sin(2 * d$x), 0.3 * d$x), 0.5)
  d$y3 <- qgamma(U[, 3], shape = 2.5, scale = exp(0.1 + d$w * cos(1.5 * d$x)) / 2.5)
  fit <- tricopula(list(y1 ~ te(x, z, k = 5), y2 ~ g + s(x, by = g), y3 ~ s(x, by = w)), d)
  expect_equal(fit$convergence, 0L)
  expect_true(fit$sp_converged)
  expect_equal(names(fit$sp), c("margin1:te(x,z).1", "margin1:te(x,z).2",
                                "margin2:s(x):ga", "margin2:s(x):gb", "margin3:s(x):w"))
  tab <- summary(fit)$smooths
  expect_equal(tab$term, c("margin1:te(x,z)", "margin2:s(x):ga", "margin2:s(x):gb",
                           "margin3:s(x):w"))
  expect_equal(tab$penalties, c(2L, 1L, 1L, 1L))
  expect_gt(tab$edf[2], 4)          # sin(2x) in group a
  expect_lt(tab$edf[3], 2.5)        # linear in group b
  expect_lt(max(abs(fit$gradient)), 1e-4)
  # penalty matrices and the two-penalty Fellner-Schall fixed point,
  # rebuilt independently
  sm <- mgcv::smoothCon(mgcv::te(x, z, k = 5), d, absorb.cons = TRUE,
                        scale.penalty = TRUE)[[1]]
  ii <- fit$smooth_blocks[[1]]$ind
  expect_equal(unname(fit$design[[1]][, -1]), unname(sm$X), tolerance = 1e-10)
  b <- coef(fit)[ii]
  St <- fit$sp[1] * sm$S[[1]] + fit$sp[2] * sm$S[[2]]
  e <- eigen(St, symmetric = TRUE); keep <- e$values > max(e$values) * 1e-8
  Sinv <- e$vectors[, keep] %*% (t(e$vectors[, keep]) / e$values[keep])
  Vp <- vcov(fit)[ii, ii]
  upd1 <- (sum(Sinv * sm$S[[1]]) - sum(Vp * sm$S[[1]])) / sum(b * (sm$S[[1]] %*% b))
  expect_equal(unname(upd1), 1, tolerance = 0.03)
  # the surface is recovered; agreement with mgcv
  g <- expand.grid(x = seq(-1.8, 1.8, length.out = 13), z = seq(-0.9, 0.9, length.out = 7))
  nd <- cbind(g, g = factor("a", levels = c("a", "b")), w = 0)
  pm <- predict(fit, nd, type = "mean")$y1
  expect_lt(sqrt(mean((pm - (1 + f2(g$x, g$z)))^2)), 0.06)
  gm <- mgcv::gam(y1 ~ te(x, z, k = 5), data = d, method = "REML")
  expect_lt(max(abs(pm - predict(gm, g))), 0.15)
  expect_equal(predict(fit, d[1:5, ], type = "mean"), predict(fit, type = "mean")[1:5, ],
               ignore_attr = TRUE, tolerance = 1e-10)
  # tc_smooth for smooths with a by variable or two covariates needs newdata
  expect_error(tc_smooth(fit, 1), "Supply newdata")
  ca <- tc_smooth(fit, "margin2:s(x):ga",
                  newdata = data.frame(x = seq(-1.5, 1.5, length.out = 15),
                                       g = factor("a", levels = c("a", "b"))))
  tr <- sin(2 * ca$x)
  expect_gt(cor(ca$fit, tr), 0.98)
  cw <- tc_smooth(fit, 4, newdata = data.frame(x = c(-1, 0, 1), w = 0))
  expect_equal(cw$fit, c(0, 0, 0))     # numeric by = 0 switches the smooth off
  s1 <- tc_smooth(fit, 1, newdata = g[1:4, ])
  expect_named(s1, c("x", "z", "fit", "se", "lower", "upper"))
  expect_equal(dim(simulate(fit, seed = 1)), c(n, 3L))
  bt <- tc_bootstrap(fit, B = 2, type = "nonparametric", seed = 8)
  expect_equal(bt$failures, 0L)
  # With independent margins the selection must agree with mgcv's REML for
  # the margin alone (under dependence the joint model legitimately selects
  # differently, because the other responses carry information).
  set.seed(177)
  di <- data.frame(x = runif(n, -2, 2), z = runif(n, -1, 1))
  di$y1 <- rnorm(n, 1 + f2(di$x, di$z), 0.5)
  di$y2 <- rnorm(n, 0.3 * di$x); di$y3 <- rnorm(n, 0.2 * di$x)
  fi <- tricopula(list(y1 ~ te(x, z, k = 5), y2 ~ x, y3 ~ x), di,
                  margins = rep("gaussian", 3), hessian = FALSE)
  gi <- mgcv::gam(y1 ~ te(x, z, k = 5), data = di, method = "REML")
  s2 <- exp(2 * coef(fi)[["margin1_log_scale"]])
  expect_equal(summary(fi)$smooths$edf, sum(gi$edf[-1]), tolerance = 0.05)
  expect_equal(unname(fi$sp[1] * s2), unname(gi$sp[1]), tolerance = 0.15)
  # REML score of mgcv at TriCopula's smoothing parameters: at its optimum
  gs <- mgcv::gam(y1 ~ te(x, z, k = 5), data = di, method = "REML", sp = fi$sp * s2)
  expect_lt(gs$gcv.ubre - gi$gcv.ubre, 0.05)
  expect_lt(max(abs(fitted(gi) - predict(fi, type = "mean")$y1)), 0.06)
  # an unpenalized smooth (fx = TRUE) has no smoothing parameter
  fx <- tricopula(list(y1 ~ s(x, k = 5, fx = TRUE), y2 ~ x, y3 ~ x), d, hessian = FALSE)
  expect_null(fx$sp); expect_null(fx$penalties)
  expect_equal(length(fx$smooth_blocks), 1L)
  expect_equal(summary(fx)$smooths$penalties, 0L)
})
