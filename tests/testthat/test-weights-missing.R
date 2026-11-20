# Prior weights and missing responses (0.13.0)

test_that("weights: defaults unchanged, integer weights equal replicated rows", {
  skip_on_cran()
  d <- sim_mixed(250, seed = 181)
  f0 <- fit_mixed(d)
  f1 <- fit_mixed(d, weights = rep(1, nrow(d)), missing = "fail")
  expect_null(f0$weights); expect_null(f0$observed)
  expect_equal(coef(f1), coef(f0), tolerance = 1e-10)
  expect_equal(f1$logLik, f0$logLik, tolerance = 1e-10)
  expect_equal(vcov(f1), vcov(f0), tolerance = 1e-8)
  set.seed(182)
  w <- sample(0:3, nrow(d), replace = TRUE)
  fw <- fit_mixed(d, weights = w)
  fr <- fit_mixed(d[rep(seq_len(nrow(d)), w), ])
  expect_equal(unname(coef(fw)), unname(coef(fr)), tolerance = 1e-5)
  expect_equal(fw$logLik, fr$logLik, tolerance = 1e-7)
  expect_equal(unname(vcov(fw)), unname(vcov(fr)), tolerance = 1e-4)
  expect_equal(unname(fw$ifm$coefficients), unname(fr$ifm$coefficients), tolerance = 1e-5)
  expect_equal(nobs(fw), nrow(d))
  expect_lt(max(abs(fw$gradient)), 1e-3)
  # log-likelihood and score at a fixed point: weighted sums of the rows
  par <- true_par()
  ll <- vapply(seq_len(nrow(d)), function(i) {
    tc_loglik(fit_mixed(d[i, ], start = par, control = list(maxit = 0), polish = FALSE,
                        hessian = FALSE), par)
  }, numeric(1))
  expect_equal(tc_loglik(fw, par), sum(w * ll), tolerance = 1e-9)
  expect_equal(unname(tc_gradient(fw, par)),
               numDeriv::grad(function(p) tc_loglik(fw, p), par), tolerance = 1e-5)
  # sandwich: rows are independent units, sum of w^2 s s'
  V <- vcov(fw); S <- fw$scores
  expect_equal(vcov(fw, type = "sandwich"), V %*% crossprod(S) %*% V, tolerance = 1e-10)
  # non-integer weights, a PCC, post-estimation and bootstrap run
  wd <- runif(nrow(d), 0.5, 2)
  fp <- fit_mixed(d, "pcc", families = c("clayton", "frank", "gaussian"), weights = wd)
  expect_equal(fp$convergence, 0L)
  expect_equal(tc_loglik(fp), sum(wd * tc_joint_logscore(fp)), tolerance = 1e-9)
  bt <- tc_bootstrap(fp, B = 3, type = "nonparametric", seed = 1)
  expect_equal(bt$failures, 0L)
  expect_error(fit_mixed(d, weights = c(1, 2)), "weights must be")
  expect_error(fit_mixed(d, weights = rep(-1, nrow(d))), "weights must be")
})

test_that("missing responses: copula terms match mvtnorm and VineCopula", {
  skip_on_cran()
  skip_if_not_installed("mvtnorm")
  set.seed(183)
  n <- 50
  U <- matrix(runif(3 * n, 0.03, 0.97), n, 3)
  U[1:9, 1] <- NA; U[10:18, 2] <- NA; U[19:27, 3] <- NA; U[28:31, 1:2] <- NA; U[32:33, 2:3] <- NA
  W <- rep(list(matrix(1, n, 1, dimnames = list(NULL, "(Intercept)"))), 3)
  th <- TriCopula:::.tc_dep_inverse(R_true)
  for (cop in c("gaussian", "student")) {
    tt <- if (cop == "student") c(th, log(5 - 2.05)) else th
    r <- TriCopula:::.tc_cop_eval(tt, U, W, cop, deriv = TRUE)
    ref <- vapply(seq_len(n), function(i) {
      o <- which(!is.na(U[i, ]))
      if (length(o) < 2) return(0)
      if (cop == "gaussian") {
        z <- qnorm(U[i, o])
        mvtnorm::dmvnorm(z, sigma = R_true[o, o], log = TRUE) - sum(dnorm(z, log = TRUE))
      } else {
        z <- qt(U[i, o], 5)
        mvtnorm::dmvt(z, sigma = R_true[o, o], df = 5, log = TRUE) - sum(dt(z, 5, log = TRUE))
      }
    }, numeric(1))
    expect_equal(r$lc, ref, tolerance = 1e-10, info = cop)
    g <- numDeriv::grad(function(x) sum(TriCopula:::.tc_cop_eval(x, U, W, cop)$lc), tt)
    expect_equal(colSums(r$score), g, tolerance = 1e-6, info = cop)
    gu <- numDeriv::grad(function(x) {
      V <- U; V[!is.na(U)] <- x
      sum(TriCopula:::.tc_cop_eval(tt, V, W, cop)$lc)
    }, U[!is.na(U)])
    expect_equal(r$du[!is.na(U)], gu, tolerance = 1e-6, info = cop)
    expect_true(all(r$du[is.na(U)] == 0))
  }
  # a Gaussian vine with any root equals the Gaussian copula for every
  # pattern (root missing: the quadrature against the closed form)
  ref <- TriCopula:::.tc_cop_eval(th, U, W, "gaussian")$lc
  for (root in 1:3) {
    o <- TriCopula:::.tc_pcc_order(root)
    cop <- structure("pcc", families = rep("gaussian", 3), root = root)
    r <- TriCopula:::.tc_cop_eval(TriCopula:::.tc_dep_inverse(R_true[o, o]), U, W, cop)
    expect_equal(r$lc, ref, tolerance = 1e-9, info = root)
    expect_true(all(is.na(r$rosenblatt[rowSums(is.na(U)) > 0, ])))
  }
  skip_if_not_installed("VineCopula")
  fam <- c("clayton", "frank", "student")
  thp <- c(log(1.5), 3, atanh(0.4), log(6 - 2.05))
  for (root in 1:3) {
    cop <- structure("pcc", families = fam, root = root)
    RVM <- VineCopula::C2RVine(TriCopula:::.tc_pcc_order(root), family = c(3, 5, 2),
                               par = c(1.5, 3, 0.4), par2 = c(0, 0, 6))
    r <- TriCopula:::.tc_cop_eval(thp, U, W, cop, deriv = TRUE)
    ref <- vapply(seq_len(n), function(i) {
      ob <- !is.na(U[i, ])
      if (sum(ob) < 2) return(0)
      if (all(ob)) return(log(VineCopula::RVinePDF(U[i, , drop = FALSE], RVM)))
      mis <- which(!ob)
      log(integrate(function(s) vapply(s, function(v) {
        u <- U[i, ]; u[mis] <- v
        VineCopula::RVinePDF(matrix(u, 1), RVM)
      }, numeric(1)), 0, 1, rel.tol = 1e-10)$value)
    }, numeric(1))
    expect_equal(r$lc, ref, tolerance = 1e-7, info = root)
    g <- numDeriv::grad(function(x) sum(TriCopula:::.tc_cop_eval(x, U, W, cop)$lc), thp)
    expect_equal(colSums(r$score), g, tolerance = 1e-5, info = root)
    gu <- numDeriv::grad(function(x) {
      V <- U; V[!is.na(U)] <- x
      sum(TriCopula:::.tc_cop_eval(thp, V, W, cop)$lc)
    }, U[!is.na(U)])
    expect_equal(r$du[!is.na(U)], gu, tolerance = 1e-5, info = root)
  }
})

hole <- function(d, seed, p = 0.25) {
  set.seed(seed)
  for (v in c("y1", "y2", "y3")) d[[v]][runif(nrow(d)) < p] <- NA
  d
}

test_that("missing = fail / complete / likelihood: data handling and likelihood", {
  skip_on_cran()
  d <- sim_mixed(400, seed = 184)
  f0 <- fit_mixed(d, hessian = FALSE)
  fl <- fit_mixed(d, hessian = FALSE, missing = "likelihood")   # nothing missing
  expect_identical(coef(fl), coef(f0))
  expect_null(fl$observed)
  dm <- hole(d, 185)
  dm$x[5] <- NA
  expect_error(fit_mixed(dm), "missing values")
  fc <- fit_mixed(dm, hessian = FALSE, missing = "complete")
  cc <- complete.cases(dm)
  expect_equal(fc$n, sum(cc)); expect_equal(fc$rows_dropped, sum(!cc))
  expect_equal(coef(fc), coef(fit_mixed(dm[cc, ], hessian = FALSE)), tolerance = 1e-10)
  fm <- fit_mixed(dm, missing = "likelihood")
  all_na <- is.na(dm$y1) & is.na(dm$y2) & is.na(dm$y3)
  expect_equal(fm$n, sum(!is.na(dm$x) & !all_na))
  expect_equal(fm$rows_dropped, sum(is.na(dm$x) | all_na))
  expect_equal(dim(fm$observed), c(fm$n, 3L))
  expect_equal(colSums(!fm$observed), colSums(is.na(fm$data[c("y1", "y2", "y3")])),
               ignore_attr = TRUE)
  expect_true(all(is.na(fm$uniforms[!fm$observed])))
  expect_equal(fm$convergence, 0L)
  expect_false(is.null(fm$vcov))
  # the observed-data log-likelihood, rebuilt from stats and mvtnorm
  skip_if_not_installed("mvtnorm")
  par <- true_par()
  dd <- fm$data
  lm <- cbind(dnorm(dd$y1, 1 + 0.5 * dd$x, 1.2, log = TRUE),
              dlnorm(dd$y2, 0.2 + 0.3 * dd$x, 0.5, log = TRUE),
              dgamma(dd$y3, shape = 2.5, scale = exp(0.1 + 0.2 * dd$x) / 2.5, log = TRUE))
  z <- qnorm(cbind(pnorm(dd$y1, 1 + 0.5 * dd$x, 1.2), plnorm(dd$y2, 0.2 + 0.3 * dd$x, 0.5),
                   pgamma(dd$y3, shape = 2.5, scale = exp(0.1 + 0.2 * dd$x) / 2.5)))
  ref <- vapply(seq_len(nrow(dd)), function(i) {
    o <- which(!is.na(z[i, ]))
    lc <- if (length(o) < 2) 0 else {
      mvtnorm::dmvnorm(z[i, o], sigma = R_true[o, o], log = TRUE) - sum(dnorm(z[i, o], log = TRUE))
    }
    sum(lm[i, o]) + lc
  }, numeric(1))
  expect_equal(tc_loglik(fm, par), sum(ref), tolerance = 1e-9)
  expect_equal(unname(tc_gradient(fm, par)),
               numDeriv::grad(function(p) tc_loglik(fm, p), par), tolerance = 1e-5)
  # more information than complete cases: smaller standard errors
  fcc <- fit_mixed(dm, missing = "complete")
  expect_true(all(sqrt(diag(vcov(fm))) < sqrt(diag(vcov(fcc)))))
  z0 <- (coef(fm) - par) / sqrt(diag(vcov(fm)))
  expect_lt(max(abs(z0)), 4)
  # post-estimation with missing responses
  expect_s3_class(tc_diagnostics(fm), "tc_diagnostics")
  expect_equal(nrow(predict(fm, type = "mean")), fm$n)
  expect_true(all(is.finite(tc_dependence(fm)$std.error)))
  bp <- tc_bootstrap(fm, B = 3, seed = 2)
  bn <- tc_bootstrap(fm, B = 3, type = "nonparametric", seed = 3)
  expect_equal(bp$failures + bn$failures, 0L)
  # weights and missing responses together
  wv <- rep(c(1, 2), length.out = nrow(dm))
  fwm <- fit_mixed(dm, missing = "likelihood", weights = wv, hessian = FALSE)
  expect_equal(length(fwm$weights), fwm$n)
  expect_equal(tc_loglik(fwm), sum(fwm$weights * tc_joint_logscore(fwm)), tolerance = 1e-9)
})

test_that("missing responses in Student-t and PCC fits, including a missing root", {
  skip_on_cran()
  fam <- c("clayton", "gumbel180", "frank")
  set.seed(186)
  n <- 700
  d <- data.frame(x = rnorm(n))
  U <- tc_rpcc(n, fam, c(2, 1.6, 3))
  d$y1 <- qnorm(U[, 1], 1 + 0.5 * d$x, 1.2)
  d$y2 <- qlnorm(U[, 2], 0.2 + 0.3 * d$x, 0.5)
  d$y3 <- qgamma(U[, 3], shape = 2.5, scale = exp(0.1 + 0.2 * d$x) / 2.5)
  dm <- hole(d, 187, 0.2)
  fit <- fit_mixed(dm, "pcc", families = fam, missing = "likelihood")
  expect_equal(fit$convergence, 0L)
  expect_false(is.null(fit$vcov))
  truth <- c(1, 0.5, log(1.2), 0.2, 0.3, log(0.5), 0.1, 0.2, log(2.5), log(2), log(0.6), 3)
  z <- (coef(fit) - truth) / sqrt(diag(vcov(fit)))
  expect_lt(max(abs(z)), 4)
  expect_lt(max(abs(fit$gradient)), 1e-2)
  expect_equal(unname(tc_gradient(fit, truth)),
               numDeriv::grad(function(p) tc_loglik(fit, p), truth), tolerance = 1e-4)
  expect_gt(sum(!fit$observed[, 1] & fit$observed[, 2] & fit$observed[, 3]), 20)
  # another root moves the integral to another pattern; the Gaussian vine
  # has the same likelihood whatever the root
  g1 <- fit_mixed(dm, "pcc", families = rep("gaussian", 3), missing = "likelihood",
                  hessian = FALSE)
  g3 <- fit_mixed(dm, "pcc", families = rep("gaussian", 3), missing = "likelihood",
                  hessian = FALSE, root = 3)
  gg <- fit_mixed(dm, missing = "likelihood", hessian = FALSE)
  expect_equal(g1$logLik, gg$logLik, tolerance = 1e-5)
  expect_equal(g3$logLik, gg$logLik, tolerance = 1e-5)
  # Student-t copula with a df predictor and missing responses
  ds <- hole(sim_mixed(600, "student", seed = 188), 189, 0.2)
  fs <- fit_mixed(ds, "student", missing = "likelihood")
  expect_equal(fs$convergence, 0L)
  zs <- (coef(fs) - true_par("student")) / sqrt(diag(vcov(fs)))
  expect_lt(max(abs(zs)), 4)
  expect_equal(unname(tc_gradient(fs, true_par("student"))),
               numDeriv::grad(function(p) tc_loglik(fs, p), true_par("student")),
               tolerance = 1e-4)
  sel <- tc_select(list(y1 ~ x, y2 ~ x, y3 ~ x), dm, missing = "likelihood", hessian = FALSE)
  expect_equal(sel$families[1], "clayton")
  expect_equal(sel$fit$missing, "likelihood")
})

test_that("weights_type: sandwich is the default covariance for sampling weights", {
  d <- sim_mixed(250, seed = 190)
  f0 <- fit_mixed(d)
  expect_null(f0$weights_type)
  expect_identical(vcov(f0), f0$vcov)                       # unweighted: as before
  set.seed(191)
  wi <- sample(1:3, nrow(d), replace = TRUE)
  fi <- fit_mixed(d, weights = wi)                          # whole numbers: frequency
  expect_equal(fi$weights_type, "frequency")
  expect_identical(vcov(fi), fi$vcov)
  expect_equal(summary(fi)$vcov_type, "hessian")
  ws <- runif(nrow(d), 0.5, 2)
  fs <- fit_mixed(d, weights = ws)                          # not whole numbers: sampling
  expect_equal(fs$weights_type, "sampling")
  V <- fs$vcov
  expect_equal(vcov(fs), V %*% crossprod(fs$scores) %*% V, tolerance = 1e-10)
  expect_identical(vcov(fs, type = "hessian"), fs$vcov)
  s <- summary(fs)
  expect_equal(s$vcov_type, "sandwich")
  expect_equal(unname(s$coefficients[, "Std. Error"]), unname(sqrt(diag(vcov(fs, type = "sandwich")))))
  expect_output(print(s), "sampling")
  ff <- fit_mixed(d, weights = ws, weights_type = "frequency")
  expect_equal(coef(ff), coef(fs))
  expect_identical(vcov(ff), ff$vcov)
  expect_equal(confint(fs), confint(ff, type = "sandwich"))
  expect_equal(tc_dependence(fs)$std.error, tc_dependence(ff, type = "sandwich")$std.error)
  # whole-number sampling weights must be declared
  fd <- fit_mixed(d, weights = wi, weights_type = "sampling")
  expect_equal(summary(fd)$vcov_type, "sandwich")
  expect_error(fit_mixed(d, weights = ws, weights_type = "other"))
})
