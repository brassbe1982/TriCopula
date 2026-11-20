test_that("prediction, simulation and dependence summaries are consistent", {
  d <- sim_mixed(400, "student", seed = 51)
  d$g <- factor(rep(c("a", "b"), 200))
  fit <- tricopula(list(y1 ~ x + g, y2 ~ x, y3 ~ x), d,
                   margins = c("gaussian", "lognormal", "gamma"), copula = "student")
  nd <- data.frame(x = c(-1, 0, 1), g = factor("b", levels = c("a", "b")))
  pr <- predict(fit, nd)
  expect_named(pr, paste0("margin", 1:3))
  expect_equal(nrow(pr$margin1), 3L)
  b <- coef(fit)
  expect_equal(pr$margin1$location,
               unname(b[1] + b[2] * nd$x + b[3]), tolerance = 1e-12)
  expect_equal(predict(fit)$margin3$mean,
               unname(fit$fitted_margins[[3]]$pars$mean))
  mu <- predict(fit, nd, type = "mean")
  expect_named(mu, c("y1", "y2", "y3"))
  expect_equal(mu$y2, exp(pr$margin2$meanlog + pr$margin2$sdlog^2 / 2))

  s1 <- simulate(fit, seed = 9)
  s2 <- simulate(fit, seed = 9)
  expect_equal(s1, s2)
  expect_equal(dim(s1), c(400L, 3L))
  expect_true(all(s1$y2 > 0 & s1$y3 > 0))
  expect_length(simulate(fit, nsim = 2, newdata = nd), 2L)

  tau <- tc_tau(fit)
  expect_equal(tau[1, 2], 2 / pi * asin(fit$correlation[1, 2]))
  lam <- tc_taildep(fit)
  expect_true(all(lam > 0 & lam <= 1))
  dep <- tc_dependence(fit)
  expect_true(all(c("rho12", "rho13", "rho23", "rho23|1", "df") %in% rownames(dep)))
  expect_true(all(dep$lower < dep$estimate & dep$estimate < dep$upper))
  expect_equal(dep["rho23", "estimate"], fit$correlation[2, 3])
  expect_equal(nrow(confint(fit)), length(b))

  fg <- tricopula(list(y1 ~ x + g, y2 ~ x, y3 ~ x), d,
                  margins = c("gaussian", "lognormal", "gamma"))
  expect_true(all(tc_taildep(fg)[upper.tri(diag(3))] == 0))
  cmp <- tc_compare(fg, fit)
  expect_equal(cmp$npar, c(13L, 14L))
  expect_lt(cmp$AIC[2], cmp$AIC[1])  # data are Student-t with df = 5
})

test_that("joint probabilities are correct", {
  skip_if_not_installed("mvtnorm")
  d <- sim_mixed(300, seed = 61)
  fit <- fit_mixed(d)
  nd <- data.frame(x = c(-0.5, 0.5))
  th <- c(1, 1.2, 1)
  jp <- tc_jointprob(fit, th, nd)
  expect_length(jp, 2L)
  expect_true(all(jp > 0 & jp < 1))
  # Monte Carlo check against simulation from the fitted model
  set.seed(62)
  sim <- simulate(fit, newdata = nd[rep(1, 2e5), , drop = FALSE])
  mc <- mean(sim$y1 <= th[1] & sim$y2 <= th[2] & sim$y3 <= th[3])
  expect_equal(jp[1], mc, tolerance = 0.02)
  expect_length(tc_jointprob(fit, th), 300L)
  expect_error(tc_jointprob(fit, 1:2), "length 3")

  # Student-t copula CDF: mixture integration vs mvtnorm at integer df
  u <- c(0.3, 0.6, 0.45)
  ref <- as.numeric(mvtnorm::pmvt(upper = qt(u, 5), corr = R_true, df = 5,
                                  abseps = 1e-6))
  expect_equal(TriCopula:::.tc_copula_cdf(u, R_true, 5, "student"), ref,
               tolerance = 1e-3)
  # independence: C(u) = prod(u)
  expect_equal(TriCopula:::.tc_copula_cdf(u, diag(3), Inf, "gaussian"), prod(u),
               tolerance = 1e-6)
})
