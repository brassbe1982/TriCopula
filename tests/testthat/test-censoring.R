# Truncation and censoring (0.16.0)

fm3 <- list(y1 ~ x, y2 ~ x, y3 ~ x)
mg3 <- c("gaussian", "lognormal", "gamma")

# censor a complete data set: y1 right-censored at a fixed point, y2 left-
# and interval-censored, y3 right-censored at random times
censor_it <- function(d, seed = 1, p3 = 0.35, c0 = 0) {
  set.seed(seed)
  n <- nrow(d)
  d$s1 <- as.integer(d$y1 <= 1.8); d$y1 <- pmin(d$y1, 1.8)
  d$s2 <- 1L; d$y2u <- NA_real_
  lf <- d$y2 < 0.8
  d$s2[lf] <- 2L; d$y2[lf] <- 0.8
  iv <- !lf & runif(n) < 0.3
  d$s2[iv] <- 3L
  d$y2u[iv] <- ceiling(d$y2[iv] * 2) / 2 + 0.01
  d$y2[iv] <- floor(d$y2[iv] * 2) / 2 + 0.01
  cc <- c0 + rexp(n, p3)
  d$s3 <- as.integer(d$y3 <= cc); d$y3 <- pmin(d$y3, cc)
  d
}
cens3 <- list(c(status = "s1"), c(status = "s2", upper = "y2u"), c(status = "s3"))

# limits on the copula scale and log-densities of the three margins at `par`
# (first nine elements), by base R functions
ref_margins <- function(d, par) {
  F1 <- function(y) pnorm(y, par[1] + par[2] * d$x, exp(par[3]))
  F2 <- function(y) plnorm(y, par[4] + par[5] * d$x, exp(par[6]))
  sh <- exp(par[9]); mu3 <- exp(par[7] + par[8] * d$x)
  F3 <- function(y) pgamma(y, shape = sh, scale = mu3 / sh)
  lo <- cbind(F1(d$y1), ifelse(d$s2 == 2L, 0, F2(d$y2)), F3(d$y3))
  up <- cbind(ifelse(d$s1 == 0L, 1, lo[, 1]),
              ifelse(d$s2 == 3L, F2(ifelse(d$s2 == 3L, d$y2u, 1)), F2(d$y2)),
              ifelse(d$s3 == 0L, 1, lo[, 3]))
  lf <- cbind(dnorm(d$y1, par[1] + par[2] * d$x, exp(par[3]), log = TRUE),
              dlnorm(d$y2, par[4] + par[5] * d$x, exp(par[6]), log = TRUE),
              dgamma(d$y3, shape = sh, scale = mu3 / sh, log = TRUE))
  cen <- cbind(d$s1 != 1L, d$s2 != 1L, d$s3 != 1L)
  lf[cen] <- 0
  list(lo = lo, up = up, cen = cen, lf = lf)
}

test_that("defaults unchanged: no truncation or censoring is the earlier model", {
  d <- sim_mixed(200, seed = 301)
  f0 <- fit_mixed(d)
  d$one <- 1L
  f1 <- fit_mixed(d, truncation = NULL, censoring = NULL)
  f2 <- fit_mixed(d, truncation = list(c(-Inf, Inf), NULL, list(lower = -Inf)),
                  censoring = list(c(status = "one"), NULL, c(status = "one")))
  # a bound at the edge of the support changes nothing but is recorded
  f3 <- fit_mixed(d, truncation = list(NULL, c(0, Inf), NULL))
  expect_equal(coef(f3), coef(f0), tolerance = 1e-8)
  expect_equal(f3$logLik, f0$logLik, tolerance = 1e-10)
  expect_null(f0$truncation); expect_null(f0$censored); expect_null(attr(f0$response, "ct"))
  for (f in list(f1, f2)) {
    expect_identical(coef(f), coef(f0))
    expect_identical(f$logLik, f0$logLik)
    expect_identical(vcov(f), vcov(f0))
    expect_null(f$truncation); expect_null(f$censored)
  }
})

test_that("truncated margins equal the formulas; type joint divides the untruncated density", {
  d <- sim_mixed(3000, seed = 302)
  d <- d[d$y2 > 1 & d$y2 < 4 & d$y1 > 0, ][1:300, ]
  par <- true_par() + 0.05
  tr <- list(c(0, Inf), c(1, 4), NULL)
  ft <- fit_mixed(d, truncation = tr, start = par, control = list(maxit = 0), polish = FALSE,
                  hessian = FALSE)
  m1 <- par[1] + par[2] * d$x; s1 <- exp(par[3])
  D1 <- 1 - pnorm(0, m1, s1)
  m2 <- par[4] + par[5] * d$x; s2 <- exp(par[6])
  D2 <- plnorm(4, m2, s2) - plnorm(1, m2, s2)
  u1 <- (pnorm(d$y1, m1, s1) - pnorm(0, m1, s1)) / D1
  u2 <- (plnorm(d$y2, m2, s2) - plnorm(1, m2, s2)) / D2
  expect_equal(unname(ft$uniforms[, 1]), u1, tolerance = 1e-10)
  expect_equal(unname(ft$uniforms[, 2]), u2, tolerance = 1e-10)
  expect_equal(unname(ft$fitted_margins[[1]]$logf), dnorm(d$y1, m1, s1, log = TRUE) - log(D1),
               tolerance = 1e-10)
  expect_equal(unname(ft$fitted_margins[[2]]$logf), dlnorm(d$y2, m2, s2, log = TRUE) - log(D2),
               tolerance = 1e-10)
  # the joint log-likelihood: truncated densities plus the copula at the
  # rescaled PITs
  f0 <- fit_mixed(d, start = par, control = list(maxit = 0), polish = FALSE, hessian = FALSE)
  W <- f0$dep_design
  lc <- TriCopula:::.tc_cop_eval(par[10:12], cbind(u1, u2, f0$uniforms[, 3]), W, "gaussian")$lc
  expect_equal(tc_loglik(ft, par),
               sum(ft$fitted_margins[[1]]$logf + ft$fitted_margins[[2]]$logf +
                     f0$fitted_margins[[3]]$logf + lc), tolerance = 1e-9)
  expect_equal(unname(tc_gradient(ft, par)),
               numDeriv::grad(function(p) tc_loglik(ft, p), par), tolerance = 1e-5)
  # type "joint": untruncated log-likelihood minus log P(inside)
  fj <- fit_mixed(d, truncation = list(NULL, c(1, 4), NULL), truncation_type = "joint",
                  start = par, control = list(maxit = 0), polish = FALSE, hessian = FALSE)
  expect_equal(tc_loglik(fj, par), tc_loglik(f0, par) - sum(log(D2)), tolerance = 1e-9)
  expect_equal(fj$uniforms, f0$uniforms, tolerance = 1e-12)
  expect_equal(unname(tc_gradient(fj, par)),
               numDeriv::grad(function(p) tc_loglik(fj, p), par), tolerance = 1e-5)
  # bounds per observation from columns of the data
  d$lo <- 1; d$hi <- 4
  fc <- fit_mixed(d, truncation = list(c(0, Inf), list(lower = "lo", upper = "hi"), NULL),
                  start = par, control = list(maxit = 0), polish = FALSE, hessian = FALSE)
  expect_equal(tc_loglik(fc, par), tc_loglik(ft, par), tolerance = 1e-12)
  expect_error(fit_mixed(d, truncation = list(NULL, c(2, 4), NULL)), "outside its truncation")
  expect_error(fit_mixed(d, truncation = list(c(0, Inf), c(1, 4), NULL),
                         truncation_type = "joint"), "one truncated response")
  expect_error(fit_mixed(d, truncation = list(c(1, 0), NULL, NULL)), "below upper")
})

test_that("censored log-likelihood, Gaussian copula: closed forms and mvtnorm", {
  skip_if_not_installed("mvtnorm")
  d <- censor_it(sim_mixed(60, seed = 303), seed = 3)
  par <- true_par() + 0.03
  f <- fit_mixed(d, censoring = cens3, start = par, control = list(maxit = 0), polish = FALSE,
                 hessian = FALSE)
  expect_true(all(tabulate(rowSums(f$censored > 0L), 3L) > 0L))   # 1, 2 and 3 censored occur
  rm <- ref_margins(d, par)
  dep <- TriCopula:::.tc_dep_vec(par[10], par[11], par[12])
  R <- TriCopula:::.tc_R(dep$r12, dep$r13, dep$r23)
  ref <- vapply(seq_len(nrow(d)), function(i) {
    K <- which(rm$cen[i, ]); E <- which(!rm$cen[i, ])
    zl <- qnorm(rm$lo[i, ]); zu <- qnorm(rm$up[i, ])
    out <- sum(rm$lf[i, ])
    if (length(E) > 1L) {      # copula density of the exact coordinates
      out <- out + mvtnorm::dmvnorm(zl[E], sigma = R[E, E], log = TRUE) - sum(dnorm(zl[E], log = TRUE))
    }
    if (!length(K)) return(out)
    if (length(E)) {           # conditional normal of the censored given the exact
      A <- R[K, E, drop = FALSE] %*% solve(R[E, E])
      mu <- drop(A %*% zl[E]); S <- R[K, K] - A %*% R[E, K, drop = FALSE]
    } else {
      mu <- rep(0, 3); S <- R
    }
    p <- if (length(K) == 1L) {
      pnorm(zu[K], mu, sqrt(S)) - pnorm(zl[K], mu, sqrt(S))
    } else {
      # (Miwa warns that it replaces infinite limits by +-1000)
      suppressWarnings(mvtnorm::pmvnorm(lower = zl[K], upper = zu[K], mean = mu,
                                        sigma = (S + t(S)) / 2,
                                        algorithm = mvtnorm::Miwa(steps = 512)))[1]
    }
    out + log(p)
  }, numeric(1))
  expect_equal(unname(tc_joint_logscore(f)), ref, tolerance = 2e-6)
  expect_equal(tc_loglik(f, par), sum(ref), tolerance = 1e-6)
  expect_equal(unname(tc_gradient(f, par)),
               numDeriv::grad(function(p) tc_loglik(f, p), par), tolerance = 1e-5)
  # more nodes change the value by less than 1e-6
  f2 <- fit_mixed(d, censoring = cens3, start = par, control = list(maxit = 0), polish = FALSE,
                  hessian = FALSE, cens_nodes = c(160, 64, 32))
  expect_equal(unname(tc_joint_logscore(f2)), ref, tolerance = 1e-6)
  expect_true(all(is.na(f$uniforms[f$censored > 0L])))
})

test_that("censored log-likelihood, vine and Student-t: adaptive integration of independent densities", {
  skip_on_cran()
  skip_if_not_installed("VineCopula")
  skip_if_not_installed("mvtnorm")
  d <- censor_it(sim_mixed(30, seed = 304), seed = 4)
  fam <- c("clayton", "gumbel", "frank"); vc <- c(3, 4, 5)
  th <- c(log(1.3), log(0.6), 2.5)
  par <- c(true_par()[1:9] + 0.03, th)
  cp <- c(1.3, 1.6, 2.5)
  for (root in c(1L, 2L)) {
    f <- fit_mixed(d, "pcc", families = fam, root = root, censoring = cens3, start = par,
                   control = list(maxit = 0), polish = FALSE, hessian = FALSE)
    rm <- ref_margins(d, par)
    ord <- TriCopula:::.tc_pcc_order(root)
    dens <- function(u) {      # vine density with VineCopula, u in variable order
      v <- u[ord]
      h1 <- VineCopula::BiCopHfunc1(v[1], v[2], vc[1], cp[1])
      h2 <- VineCopula::BiCopHfunc1(v[1], v[3], vc[2], cp[2])
      VineCopula::BiCopPDF(v[1], v[2], vc[1], cp[1]) * VineCopula::BiCopPDF(v[1], v[3], vc[2], cp[2]) *
        VineCopula::BiCopPDF(h1, h2, vc[3], cp[3])
    }
    rows <- which(rowSums(rm$cen) %in% 1:2)
    rows <- c(rows[rowSums(rm$cen)[rows] == 1L][1:3],
              if (root == 2L) rows[rowSums(rm$cen)[rows] == 2L][1])
    rows <- rows[!is.na(rows)]
    ref <- vapply(rows, function(i) {
      K <- which(rm$cen[i, ])
      g1 <- function(s, fix = NULL) vapply(s, function(ss) {
        u <- rm$lo[i, ]; u[K[1]] <- ss
        if (!is.null(fix)) u[K[2]] <- fix
        dens(u)
      }, numeric(1))
      val <- if (length(K) == 1L) {
        integrate(g1, rm$lo[i, K], rm$up[i, K], rel.tol = 1e-9, subdivisions = 300L)$value
      } else {
        integrate(function(t) vapply(t, function(tt) {
          integrate(g1, rm$lo[i, K[1]], rm$up[i, K[1]], fix = tt, rel.tol = 1e-8,
                    subdivisions = 300L)$value
        }, numeric(1)), rm$lo[i, K[2]], rm$up[i, K[2]], rel.tol = 1e-7)$value
      }
      sum(rm$lf[i, ]) + log(val)
    }, numeric(1))
    expect_equal(unname(tc_joint_logscore(f)[rows]), ref, tolerance = 1e-5)
    expect_equal(unname(tc_gradient(f, par)),
                 numDeriv::grad(function(p) tc_loglik(f, p), par), tolerance = 1e-5)
  }
  # Student-t copula, one censored response
  par <- c(true_par("student", 5)) + 0.03
  f <- fit_mixed(d, "student", censoring = cens3, start = par, control = list(maxit = 0),
                 polish = FALSE, hessian = FALSE)
  rm <- ref_margins(d, par)
  dep <- TriCopula:::.tc_dep_vec(par[10], par[11], par[12])
  R <- TriCopula:::.tc_R(dep$r12, dep$r13, dep$r23)
  nu <- 2.05 + exp(par[13])
  rows <- which(rowSums(rm$cen) == 1L)[1:6]
  ref <- vapply(rows, function(i) {
    K <- which(rm$cen[i, ])
    g <- function(s) vapply(s, function(ss) {
      u <- rm$lo[i, ]; u[K] <- ss
      z <- qt(u, nu)
      exp(mvtnorm::dmvt(z, sigma = R, df = nu, log = TRUE) - sum(dt(z, nu, log = TRUE)))
    }, numeric(1))
    sum(rm$lf[i, ]) + log(integrate(g, rm$lo[i, K], rm$up[i, K], rel.tol = 1e-9)$value)
  }, numeric(1))
  expect_equal(unname(tc_joint_logscore(f)[rows]), ref, tolerance = 1e-6)
  expect_equal(unname(tc_gradient(f, par)),
               numDeriv::grad(function(p) tc_loglik(f, p), par), tolerance = 1e-4)
})

test_that("censored marginal stage equals survival::survreg", {
  skip_on_cran()
  skip_if_not_installed("survival")
  d <- sim_mixed(300, seed = 305)
  set.seed(5)
  cc <- rexp(300, 0.4)
  d$s2 <- as.integer(d$y2 <= cc); d$y2 <- pmin(d$y2, cc)
  d$s1 <- as.integer(d$y1 <= 2); d$y1 <- pmin(d$y1, 2)
  f <- fit_mixed(d, censoring = list(c(status = "s1"), c(status = "s2"), NULL), hessian = FALSE)
  s1 <- survival::survreg(survival::Surv(y1, s1) ~ x, data = d, dist = "gaussian")
  s2 <- survival::survreg(survival::Surv(y2, s2) ~ x, data = d, dist = "lognormal")
  ifm <- f$ifm$coefficients
  expect_equal(unname(ifm[1:3]), unname(c(coef(s1), log(s1$scale))), tolerance = 1e-4)
  expect_equal(unname(ifm[4:6]), unname(c(coef(s2), log(s2$scale))), tolerance = 1e-4)
  expect_equal(f$convergence, 0L)
  expect_gte(f$logLik, f$ifm$logLik - 1e-6)
  expect_lt(max(abs(f$gradient)), 1e-3)
})

test_that("censoring and truncation together: fit, score and post-estimation", {
  skip_on_cran()
  d <- sim_mixed(2500, seed = 306)
  d <- d[d$y3 > 0.5, ][1:220, ]
  d <- censor_it(d, seed = 6, p3 = 0.15, c0 = 0.5)
  tr <- list(NULL, NULL, c(0.5, Inf))
  for (tt in c("margin", "joint")) {
    f <- fit_mixed(d, censoring = cens3, truncation = tr, truncation_type = tt,
                   cens_nodes = c(32, 16, 10))
    expect_equal(f$convergence, 0L)
    expect_lt(max(abs(f$gradient)), 1e-2)
    p <- coef(f) + 0.02
    expect_equal(unname(tc_gradient(f, p)), numDeriv::grad(function(q) tc_loglik(f, q), p),
                 tolerance = 1e-5)
    expect_true(all(is.finite(sqrt(diag(vcov(f))))))
    expect_output(print(f), "Censored values per response")
    expect_output(print(summary(f)), paste0("Truncation \\(type ", tt))
    expect_equal(dim(tc_pit(f)), c(220L, 3L))
    sim <- simulate(f, seed = 1)
    expect_true(all(sim$y3 >= 0.5))
  }
  bt <- tc_bootstrap(f, B = 2, type = "nonparametric", seed = 1)
  expect_equal(nrow(bt$replicates), 2L)
  expect_error(tc_bootstrap(f, B = 2), "parametric bootstrap is not available")
  expect_error(tc_select(fm3, d, margins = mg3, censoring = cens3), "does not support")
  d$bad <- 5L
  expect_error(fit_mixed(d, censoring = list(c(status = "bad"), NULL, NULL)), "status must be")
  expect_error(fit_mixed(d, censoring = list(NULL, c(status = "s2"), NULL)), "need upper")
  dm <- d; dm$y1[1] <- NA
  expect_error(fit_mixed(dm, censoring = cens3, missing = "likelihood"), "cannot be combined")
})

test_that("truncation: means, simulation and joint probabilities", {
  skip_on_cran()
  d <- sim_mixed(4000, seed = 307)
  d <- d[d$y2 > 1, ][1:400, ]
  tr <- list(NULL, c(1, Inf), NULL)
  fmg <- fit_mixed(d, truncation = tr)
  fjt <- fit_mixed(d, truncation = tr, truncation_type = "joint")
  nd <- data.frame(x = 0.3)
  # mean of the truncated log-normal against integrate()
  b <- coef(fmg); ml <- b[4] + b[5] * 0.3; sl <- exp(b[6])
  ref <- integrate(function(y) y * dlnorm(y, ml, sl), 1, Inf, rel.tol = 1e-10)$value /
    (1 - plnorm(1, ml, sl))
  expect_equal(predict(fmg, nd, type = "mean")[[2]], ref, tolerance = 1e-6)
  expect_equal(predict(fjt, nd, type = "mean")[[2]],
               unname(exp(coef(fjt)[4] + coef(fjt)[5] * 0.3 + exp(coef(fjt)[6])^2 / 2)),
               tolerance = 1e-10)
  skip_if_not_installed("mvtnorm")
  thr <- c(1.5, 1.8, 1.2)
  big <- nd[rep(1, 40000), , drop = FALSE]
  for (f in list(fmg, fjt)) {
    sim <- simulate(f, seed = 2, newdata = big)
    expect_true(all(sim$y2 >= 1))
    emp <- mean(sim$y1 <= thr[1] & sim$y2 <= thr[2] & sim$y3 <= thr[3])
    expect_equal(tc_jointprob(f, thr, nd), emp, tolerance = 0.01)
  }
})

# 0.17.0: closed form for one censored coordinate, accuracy check ------------

test_that("cens_method = \"quadrature\" is the censored likelihood of 0.16.0; the default is \"closed\"", {
  skip_on_cran()
  fx <- readRDS(test_path("fixtures", "censoring_0.16.0.rds"))   # computed with 0.16.0
  d <- censor_it(sim_mixed(60, seed = 303), seed = 3)
  par <- true_par() + 0.03
  f <- fit_mixed(d, censoring = cens3, start = par, control = list(maxit = 0), polish = FALSE,
                 hessian = FALSE, cens_method = "quadrature")
  expect_equal(unname(f$loglik_i), fx$gauss$ll_i, tolerance = 1e-12)
  expect_equal(unname(tc_gradient(f, par)), fx$gauss$grad, tolerance = 1e-10)
  par2 <- c(true_par()[1:9] + 0.03, log(1.3), log(0.6), 2.5)
  d2 <- censor_it(sim_mixed(30, seed = 304), seed = 4)
  f2 <- fit_mixed(d2, "pcc", families = c("clayton", "gumbel", "frank"), root = 2L,
                  censoring = cens3, start = par2, control = list(maxit = 0), polish = FALSE,
                  hessian = FALSE, cens_method = "quadrature")
  expect_equal(unname(f2$loglik_i), fx$pcc$ll_i, tolerance = 1e-12)
  expect_equal(unname(tc_gradient(f2, par2)), fx$pcc$grad, tolerance = 1e-10)
  # a complete fit
  d3 <- censor_it(sim_mixed(80, seed = 308), seed = 8)
  f3 <- fit_mixed(d3, censoring = cens3, cens_nodes = c(24, 12, 8), cens_method = "quadrature")
  expect_equal(coef(f3), fx$fit$coef, tolerance = 1e-9)
  expect_equal(f3$logLik, fx$fit$logLik, tolerance = 1e-11)
  expect_equal(vcov(f3), fx$fit$vcov, tolerance = 1e-7)
  expect_identical(f3$cens_method, "quadrature")
  expect_s3_class(f3$cens_accuracy, "tc_cens_check")
  # the check does not change the estimates
  f4 <- fit_mixed(d3, censoring = cens3, cens_nodes = c(24, 12, 8), cens_method = "quadrature",
                  cens_check = FALSE)
  expect_identical(coef(f4), coef(f3))
  expect_identical(f4$logLik, f3$logLik)
  expect_identical(vcov(f4), vcov(f3))
  expect_null(f4$cens_accuracy)
  # the default (0.17.1): closed forms; the same model, so the same estimates
  # up to the quadrature error of these few nodes
  f5 <- fit_mixed(d3, censoring = cens3)
  expect_identical(f5$cens_method, "closed")
  expect_equal(coef(f5), coef(f3), tolerance = 0.02)
  expect_equal(f5$logLik, f3$logLik, tolerance = 1e-3)
})

test_that("d C / d v of the pair copulas equals VineCopula::BiCopHfunc2", {
  skip_if_not_installed("VineCopula")
  code <- c(gaussian = 1, student = 2, clayton = 3, gumbel = 4, frank = 5, joe = 6)
  set.seed(31)
  u <- runif(40, 0.02, 0.98); v <- runif(40, 0.02, 0.98)
  for (fam in TriCopula:::.tc_pair_families) {
    base <- TriCopula:::.tc_pair_base(fam); rot <- TriCopula:::.tc_pair_rot(fam)
    eta <- if (base == "frank") 3 else 0.4
    tau <- if (base == "student") log(5 - 2.05)
    th <- TriCopula:::.tc_pair_natural(eta, fam)
    vc <- unname(code[base]) + c(`0` = 0, `180` = 10, `90` = 20, `270` = 30)[[as.character(rot)]]
    ref <- VineCopula::BiCopHfunc2(u, v, family = vc, par = if (rot %in% c(90, 270)) -th else th,
                                   par2 = if (base == "student") 5 else 0)
    expect_equal(TriCopula:::.tc_pair_h2(u, v, eta, fam, tau), ref, tolerance = 1e-9, info = fam)
  }
})

test_that("cens_method = \"closed\", Gaussian and Student-t: closed forms, mvtnorm, integrate()", {
  skip_on_cran()
  skip_if_not_installed("mvtnorm")
  d <- censor_it(sim_mixed(60, seed = 303), seed = 3)
  par <- true_par() + 0.03
  f <- fit_mixed(d, censoring = cens3, start = par, control = list(maxit = 0), polish = FALSE,
                 hessian = FALSE, cens_method = "closed")
  rm <- ref_margins(d, par)
  dep <- TriCopula:::.tc_dep_vec(par[10], par[11], par[12])
  R <- TriCopula:::.tc_R(dep$r12, dep$r13, dep$r23)
  ref <- vapply(seq_len(nrow(d)), function(i) {
    K <- which(rm$cen[i, ]); E <- which(!rm$cen[i, ])
    zl <- qnorm(rm$lo[i, ]); zu <- qnorm(rm$up[i, ])
    out <- sum(rm$lf[i, ])
    if (length(E) > 1L) {
      out <- out + mvtnorm::dmvnorm(zl[E], sigma = R[E, E], log = TRUE) - sum(dnorm(zl[E], log = TRUE))
    }
    if (!length(K)) return(out)
    if (length(E)) {
      A <- R[K, E, drop = FALSE] %*% solve(R[E, E])
      mu <- drop(A %*% zl[E]); S <- R[K, K] - A %*% R[E, K, drop = FALSE]
    } else {
      mu <- rep(0, 3); S <- R
    }
    p <- if (length(K) == 1L) {
      pnorm(zu[K], mu, sqrt(S)) - pnorm(zl[K], mu, sqrt(S))
    } else {
      suppressWarnings(mvtnorm::pmvnorm(lower = zl[K], upper = zu[K], mean = mu,
                                        sigma = (S + t(S)) / 2,
                                        algorithm = mvtnorm::Miwa(steps = 512)))[1]
    }
    out + log(p)
  }, numeric(1))
  one <- rowSums(rm$cen) == 1L
  expect_gt(sum(one), 5L)
  # one censored response: no quadrature at all
  expect_equal(unname(tc_joint_logscore(f))[one], ref[one], tolerance = 1e-11)
  # two and three: limited by the reference (Miwa)
  expect_equal(unname(tc_joint_logscore(f)), ref, tolerance = 2e-6)
  expect_equal(unname(tc_gradient(f, par)),
               numDeriv::grad(function(p) tc_loglik(f, p), par), tolerance = 1e-5)
  # Student-t copula, one censored response, against integrate() of mvtnorm::dmvt
  par <- c(true_par("student", 5)) + 0.03
  d <- censor_it(sim_mixed(30, seed = 304), seed = 4)
  f <- fit_mixed(d, "student", censoring = cens3, start = par, control = list(maxit = 0),
                 polish = FALSE, hessian = FALSE, cens_method = "closed")
  rm <- ref_margins(d, par)
  dep <- TriCopula:::.tc_dep_vec(par[10], par[11], par[12])
  R <- TriCopula:::.tc_R(dep$r12, dep$r13, dep$r23)
  nu <- 2.05 + exp(par[13])
  rows <- which(rowSums(rm$cen) == 1L)
  expect_setequal(unique(apply(rm$cen[rows, ], 1, which)), 1:3)   # each response is censored somewhere
  ref <- vapply(rows, function(i) {
    K <- which(rm$cen[i, ])
    g <- function(s) vapply(s, function(ss) {
      u <- rm$lo[i, ]; u[K] <- ss
      z <- qt(u, nu)
      exp(mvtnorm::dmvt(z, sigma = R, df = nu, log = TRUE) - sum(dt(z, nu, log = TRUE)))
    }, numeric(1))
    sum(rm$lf[i, ]) + log(integrate(g, rm$lo[i, K], rm$up[i, K], rel.tol = 1e-11)$value)
  }, numeric(1))
  expect_equal(unname(tc_joint_logscore(f)[rows]), ref, tolerance = 1e-8)
  expect_equal(unname(tc_gradient(f, par)),
               numDeriv::grad(function(p) tc_loglik(f, p), par), tolerance = 1e-4)
})

test_that("cens_method = \"closed\", vine with rotated pairs and every root: integrate() of VineCopula densities", {
  skip_on_cran()
  skip_if_not_installed("VineCopula")
  d <- censor_it(sim_mixed(40, seed = 309), seed = 9)
  fam <- c("clayton90", "gumbel", "joe270"); vc <- c(23, 4, 36)
  th <- c(log(1.3), log(0.6), log(0.5))
  cp <- c(-1.3, 1.6, -1.5)
  par <- c(true_par()[1:9] + 0.03, th)
  for (root in 1:3) {
    f <- fit_mixed(d, "pcc", families = fam, root = root, censoring = cens3, start = par,
                   control = list(maxit = 0), polish = FALSE, hessian = FALSE,
                   cens_method = "closed")
    rm <- ref_margins(d, par)
    ord <- TriCopula:::.tc_pcc_order(root)
    dens <- function(u) {
      v <- u[ord]
      h1 <- VineCopula::BiCopHfunc1(v[1], v[2], vc[1], cp[1])
      h2 <- VineCopula::BiCopHfunc1(v[1], v[3], vc[2], cp[2])
      VineCopula::BiCopPDF(v[1], v[2], vc[1], cp[1]) * VineCopula::BiCopPDF(v[1], v[3], vc[2], cp[2]) *
        VineCopula::BiCopPDF(h1, h2, vc[3], cp[3])
    }
    nc <- rowSums(rm$cen)
    # one row for each single censored response (root and both leaves), and
    # for the first root one row with two censored responses
    rows <- c(vapply(1:3, function(k) which(nc == 1L & rm$cen[, k])[1], integer(1)),
              if (root == 1L) which(nc == 2L)[1])
    expect_false(anyNA(rows))
    ref <- vapply(rows, function(i) {
      K <- which(rm$cen[i, ])
      g1 <- function(s, fix = NULL) vapply(s, function(ss) {
        u <- rm$lo[i, ]; u[K[1]] <- ss
        if (!is.null(fix)) u[K[2]] <- fix
        dens(u)
      }, numeric(1))
      val <- if (length(K) == 1L) {
        integrate(g1, rm$lo[i, K], rm$up[i, K], rel.tol = 1e-9, subdivisions = 300L)$value
      } else {
        integrate(function(t) vapply(t, function(tt) {
          integrate(g1, rm$lo[i, K[1]], rm$up[i, K[1]], fix = tt, rel.tol = 1e-7,
                    subdivisions = 300L)$value
        }, numeric(1)), rm$lo[i, K[2]], rm$up[i, K[2]], rel.tol = 1e-6)$value
      }
      sum(rm$lf[i, ]) + log(val)
    }, numeric(1))
    expect_equal(unname(tc_joint_logscore(f)[rows[1:3]]), ref[1:3], tolerance = 1e-8)
    if (root == 1L) expect_equal(unname(tc_joint_logscore(f)[rows[4]]), ref[4], tolerance = 1e-5)
    expect_equal(unname(tc_gradient(f, par)),
                 numDeriv::grad(function(p) tc_loglik(f, p), par), tolerance = 1e-5)
  }
})

test_that("cens_method = \"closed\": non-simplified vine, fits and bootstrap agree with quadrature", {
  skip_on_cran()
  d <- censor_it(sim_mixed(150, seed = 310), seed = 10)
  # non-simplified vine at fixed parameters: closed form against fine quadrature
  fam <- c("gaussian", "frank", "gumbel")
  par <- c(true_par()[1:9] + 0.03, 0.5, 2, log(0.4), 0.4, -0.3, 0.2)
  ds <- censor_it(sim_mixed(50, seed = 311), seed = 11)
  fit0 <- function(...) fit_mixed(ds, "pcc", families = fam, simplified = FALSE, k_cond = 4,
                                  censoring = cens3, start = par, sp = 1,
                                  control = list(maxit = 0), polish = FALSE, hessian = FALSE,
                                  cens_check = FALSE, ...)
  fc <- fit0(cens_method = "closed")
  # contributions at `par` (a penalized fit moves away from `start` even
  # with maxit = 0), with a given method and node numbers
  at_par <- function(method, nodes) {
    y <- fc$response
    ct <- attr(y, "ct"); ct$method <- method; ct$nodes <- as.integer(nodes)
    attr(y, "ct") <- ct
    TriCopula:::.tc_eval(par, y, fc$design, fc$margins, TriCopula:::.tc_cop(fc), fc$index,
                         W = fc$dep_design, Z = fc$scale_design, V = TriCopula:::.tc_V(fc),
                         Q = fc$shape_design)$ll_i
  }
  lcl <- at_par("closed", c(64, 32, 16))
  lq <- at_par("quadrature", c(128, 48, 16))
  nc <- rowSums(fc$censored > 0L)
  # root exact: closed form, and quadrature over at most one smooth coordinate
  leaf <- nc <= 2L & fc$censored[, 1] == 0L
  expect_gt(sum(leaf & nc > 0L), 5L)
  expect_equal(lcl[leaf], lq[leaf], tolerance = 1e-7)
  # root censored (quadrature over the root in both methods) and three
  # censored: limited by the quadrature of the reference
  expect_equal(lcl[nc <= 2L], lq[nc <= 2L], tolerance = 1e-4)
  expect_equal(lcl, lq, tolerance = 1e-2)
  expect_equal(unname(tc_gradient(fc, par)),
               numDeriv::grad(function(p) tc_loglik(fc, p), par), tolerance = 1e-5)
  # complete fits
  gq <- fit_mixed(d, censoring = cens3, cens_method = "quadrature")
  gc <- fit_mixed(d, censoring = cens3)
  expect_equal(gc$convergence, 0L)
  expect_equal(coef(gc), coef(gq), tolerance = 1e-4)
  expect_equal(gc$logLik, gq$logLik, tolerance = 1e-6)
  expect_equal(sqrt(diag(vcov(gc))), sqrt(diag(vcov(gq))), tolerance = 1e-3)
  expect_identical(gc$cens_method, "closed")
  bt <- tc_bootstrap(gc, B = 2, type = "nonparametric", seed = 1)
  expect_equal(nrow(bt$replicates), 2L)
  expect_error(fit_mixed(d, censoring = cens3, cens_method = "exact"), "should be one of")
})

test_that("tc_cens_check() measures the quadrature error and tricopula() warns", {
  skip_on_cran()
  d <- censor_it(sim_mixed(150, seed = 310), seed = 10)
  gq <- fit_mixed(d, censoring = cens3, cens_method = "quadrature")
  ck <- gq$cens_accuracy
  expect_named(ck, c("table", "loglik_change", "shift", "max_shift", "nodes", "nodes_reference",
                     "method"))
  expect_identical(ck$nodes_reference, c(128L, 64L, 32L))
  expect_equal(ck$table$rows, tabulate(rowSums(gq$censored > 0L), 3L))
  expect_lt(ck$max_shift, 0.01)
  expect_length(ck$shift, length(coef(gq)))
  expect_output(print(ck), "Largest first-order change of an estimate")
  # the first-order prediction against an actual refit with the finer nodes
  coarse <- suppressWarnings(fit_mixed(d, censoring = cens3, cens_nodes = c(8, 6, 4),
                                       cens_method = "quadrature"))
  fine <- fit_mixed(d, censoring = cens3, cens_nodes = c(16, 12, 8), cens_check = FALSE,
                    cens_method = "quadrature")
  ckc <- tc_cens_check(coarse)
  expect_identical(ckc$nodes_reference, c(16L, 12L, 8L))
  actual <- (coef(fine) - coef(coarse)) / sqrt(diag(vcov(coarse)))
  expect_equal(unname(ckc$shift), unname(actual), tolerance = 0.2)
  expect_lt(ckc$max_shift, 0.1)
  # strong dependence and few nodes: the fit warns
  set.seed(5)
  ds <- data.frame(x = rnorm(150))
  U <- tc_rcopula(150, matrix(c(1, .9, .8, .9, 1, .85, .8, .85, 1), 3))
  ds$y1 <- qnorm(U[, 1], 1 + 0.5 * ds$x, 1.2); ds$y2 <- qlnorm(U[, 2], 0.2 + 0.3 * ds$x, 0.5)
  ds$y3 <- qgamma(U[, 3], shape = 2.5, scale = exp(0.1 + 0.2 * ds$x) / 2.5)
  ds <- censor_it(ds, seed = 11)
  expect_warning(bad <- fit_mixed(ds, censoring = cens3, cens_nodes = c(8, 6, 4),
                                  cens_method = "quadrature"),
                 "not accurate enough")
  expect_gt(bad$cens_accuracy$max_shift, 0.1)
  expect_silent(good <- fit_mixed(ds, censoring = cens3, cens_nodes = c(8, 6, 4)))
  expect_lt(good$cens_accuracy$max_shift, 0.1)
  # a closed-form fit is also compared with pure quadrature
  gc <- fit_mixed(d, censoring = cens3)
  expect_lt(max(gc$cens_accuracy$table$max_diff_quadrature[1:2]), 1e-6)
  expect_lt(gc$cens_accuracy$max_shift, 1e-3)
  expect_error(tc_cens_check(fit_mixed(sim_mixed(100, seed = 1))), "no censored responses")
  expect_error(tc_cens_check(gq, nodes = c(1, 2)), "three numbers")
})
