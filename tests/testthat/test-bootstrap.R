test_that("bootstrap standard errors are reproducible and close to model-based ones", {
  skip_on_cran()
  d <- sim_mixed(300, seed = 91)
  fit <- fit_mixed(d)
  b1 <- tc_bootstrap(fit, B = 40, seed = 1)
  b2 <- tc_bootstrap(fit, B = 40, seed = 1)
  expect_equal(b1$replicates, b2$replicates)
  expect_s3_class(b1, "tc_bootstrap")
  expect_equal(b1$failures, 0L)
  expect_equal(nrow(b1$coefficients), length(coef(fit)))
  expect_true(all(c("rho12", "rho23|1") %in% rownames(b1$natural)))
  ratio <- b1$coefficients$boot.se / b1$coefficients$model.se
  expect_true(all(ratio > 0.6 & ratio < 1.6))
  expect_true(all(b1$coefficients$lower < b1$coefficients$upper))
  expect_output(print(b1), "parametric bootstrap")

  bn <- tc_bootstrap(fit, B = 20, type = "nonparametric", seed = 2)
  expect_equal(bn$type, "nonparametric")
  expect_true(all(is.finite(bn$coefficients$boot.se)))
})

test_that("bootstrap works with covariate-dependent scale and dependence, and checks inputs", {
  skip_on_cran()
  d <- sim_mixed(300, "student", seed = 92)
  fit <- fit_mixed(d, "student", scale = list(~ x, ~ 1, ~ 1), dependence = list(~ x, ~ 1, ~ 1))
  b <- tc_bootstrap(fit, B = 10, seed = 3)
  expect_equal(ncol(b$replicates), length(coef(fit)) + nrow(tc_dependence(fit)))
  d$ly <- d$y2
  f2 <- tricopula(list(y1 ~ x, log(ly) ~ x, y3 ~ x), d,
                  margins = c("gaussian", "gaussian", "gamma"), hessian = FALSE)
  expect_error(tc_bootstrap(f2, B = 5), "plain columns")
  expect_error(tc_bootstrap(1), "tricopula fit")
})
