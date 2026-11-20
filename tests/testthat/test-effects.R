test_that("average marginal effects on means have the known closed forms", {
  d <- sim_mixed(500, seed = 101)
  fit <- fit_mixed(d)
  ef <- tc_effects(fit, "x")
  b <- coef(fit); se <- sqrt(diag(vcov(fit)))
  expect_equal(ef$outcome, c("y1", "y2", "y3"))
  # Gaussian margin: effect = slope, SE = SE of the slope
  expect_equal(ef$effect[1], unname(b["margin1_x"]), tolerance = 1e-6)
  expect_equal(ef$std.error[1], unname(se["margin1_x"]), tolerance = 1e-5)
  # log-normal and Gamma with constant scale: slope times average mean
  mu <- predict(fit, type = "mean")
  expect_equal(ef$effect[2], unname(b["margin2_x"]) * mean(mu$y2), tolerance = 1e-5)
  expect_equal(ef$effect[3], unname(b["margin3_x"]) * mean(mu$y3), tolerance = 1e-5)
  expect_true(all(ef$lower < ef$effect & ef$effect < ef$upper))
  expect_true(all(ef$p.value < 0.05))

  # discrete change
  ed <- tc_effects(fit, "x", values = c(0, 1))
  nd0 <- transform(d, x = 0); nd1 <- transform(d, x = 1)
  expect_equal(ed$effect, unname(colMeans(predict(fit, nd1, type = "mean") -
                                            predict(fit, nd0, type = "mean"))))
  expect_equal(attr(ed, "kind"), "average change from 0 to 1")
  expect_equal(dim(tc_effects(fit, "x", vcov_type = "sandwich")), dim(ef))
})

test_that("effects work for factors, joint probabilities and bad inputs", {
  skip_if_not_installed("mvtnorm")
  d <- sim_mixed(400, seed = 102)
  d$g <- factor(rep(c("a", "b"), 200))
  d$y1 <- d$y1 + 0.8 * (d$g == "b")
  fit <- tricopula(list(y1 ~ x + g, y2 ~ x, y3 ~ x), d,
                   margins = c("gaussian", "lognormal", "gamma"), dependence = ~ x)
  eg <- tc_effects(fit, "g", values = c("a", "b"))
  expect_equal(eg$effect[1], unname(coef(fit)["margin1_gb"]), tolerance = 1e-8)
  expect_equal(eg$effect[2:3], c(0, 0), tolerance = 1e-10)
  expect_error(tc_effects(fit, "g"), "values")
  expect_error(tc_effects(fit, "g", values = c("a", "z")), "valid levels")
  expect_error(tc_effects(fit, "nope"), "name of one column")

  nd <- d[1:6, ]
  th <- c(1, 1.2, 1)
  ej <- tc_effects(fit, "x", newdata = nd, type = "jointprob", thresholds = th, h = 1e-3)
  ref <- mean(tc_jointprob(fit, th, transform(nd, x = x + 1e-3)) -
                tc_jointprob(fit, th, transform(nd, x = x - 1e-3))) / 2e-3
  expect_equal(ej$effect, ref, tolerance = 1e-8)
  expect_lt(ej$effect, 0)          # every mean increases with x
  expect_true(is.finite(ej$std.error) && ej$std.error > 0)
  expect_error(tc_effects(fit, "x", type = "jointprob"), "thresholds")
})
