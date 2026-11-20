test_that("Gaussian and Student copula simulation is reproducible", {
  a <- tc_simulate(100, seed = 42)
  b <- tc_simulate(100, seed = 42)
  expect_equal(a, b)
  expect_equal(dim(a), c(100L, 3L))
  s <- tc_simulate(100, copula = "student", df = 5, seed = 7)
  expect_equal(dim(s), c(100L, 3L))
  expect_error(tc_rcopula(5, diag(2)), "positive-definite")
})

test_that("copula log-densities match multivariate density identities", {
  skip_if_not_installed("mvtnorm")
  set.seed(3)
  U <- matrix(runif(60), 20, 3)
  lg <- TriCopula:::.tc_logcop(U, R_true, Inf, "gaussian")$lc
  z <- qnorm(U)
  expect_equal(lg, mvtnorm::dmvnorm(z, sigma = R_true, log = TRUE) -
                 rowSums(dnorm(z, log = TRUE)), tolerance = 1e-10)
  for (df in c(2.5, 7.3)) {
    lt <- TriCopula:::.tc_logcop(U, R_true, df, "student")$lc
    z <- qt(U, df)
    expect_equal(lt, mvtnorm::dmvt(z, sigma = R_true, df = df, log = TRUE) -
                   rowSums(dt(z, df, log = TRUE)), tolerance = 1e-10)
  }
})

test_that("Student copula density approaches the Gaussian one as df grows", {
  set.seed(4)
  U <- matrix(runif(30, 0.05, 0.95), 10, 3)
  expect_equal(TriCopula:::.tc_logcop(U, R_true, 1e7, "student")$lc,
               TriCopula:::.tc_logcop(U, R_true, Inf, "gaussian")$lc,
               tolerance = 1e-5)
})

test_that("every unconstrained parameter gives a valid correlation matrix", {
  set.seed(5)
  for (i in 1:200) {
    th <- rnorm(4, sd = 3)
    dep <- TriCopula:::.tc_dep(th, "student")
    expect_gt(min(eigen(dep$R, symmetric = TRUE)$values), 0)
    expect_gt(dep$df, 2.05)
  }
  th <- c(0.4, -0.7, 1.1, 0.3)
  dep <- TriCopula:::.tc_dep(th, "student")
  expect_equal(TriCopula:::.tc_dep_inverse(dep$R, dep$df), th, tolerance = 1e-10)
})

test_that("copula densities integrate to one (Monte Carlo)", {
  set.seed(6)
  U <- matrix(runif(3e5), ncol = 3)
  expect_equal(mean(exp(TriCopula:::.tc_logcop(U, R_true, Inf, "gaussian")$lc)),
               1, tolerance = 0.05)
  expect_equal(mean(exp(TriCopula:::.tc_logcop(U, R_true, 6, "student")$lc)),
               1, tolerance = 0.05)
})
