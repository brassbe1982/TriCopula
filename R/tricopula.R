#' TriCopula: trivariate continuous copula regression (research prototype)
#'
#' Full maximum-likelihood estimation of three continuous marginal regression
#' equations (Gaussian, log-normal or Gamma, freely mixed) linked by a
#' trivariate Gaussian or Student-t copula. All marginal and dependence
#' parameters are estimated in one optimization of the joint log-likelihood
#' \deqn{\ell = \sum_i \{\sum_{j=1}^3 \log f_j(y_{ji}) + \log c(u_{1i}, u_{2i}, u_{3i})\},
#'   \quad u_{ji} = F_j(y_{ji}).}
#'
#' This differs from the sequential workflow of `GJRM::gjrm.pcc()`, which
#' takes three previously fitted `GJRM::gamlss()` margins as given and then
#' estimates a pair-copula construction.
#'
#' @name TriCopula-package
#' @keywords internal
#' @importFrom stats simulate nobs logLik vcov confint predict
NULL

.tc_make_X <- function(formula, data, basis = NULL, na_ok = FALSE) {
  if (!inherits(formula, "formula")) stop("Each margin specification must be a formula.")
  # s() terms (0.12.0): the parametric part goes through model.frame as
  # before, the smooth bases are appended
  st <- if (.tc_has_smooth(formula)) .tc_smooth_setup(formula, data, basis)
  # na_ok: missing responses are allowed (covariates were checked before)
  mf <- stats::model.frame(if (is.null(st)) formula else st$pf, data = data,
                           na.action = if (na_ok) stats::na.pass else stats::na.fail)
  y <- stats::model.response(mf)
  if (!is.numeric(y) || is.matrix(y) || is.factor(y)) {
    stop("Each response must be a numeric continuous vector.")
  }
  tt <- stats::terms(mf)
  X <- stats::model.matrix(tt, data = mf)
  ct <- attr(X, "contrasts")
  sm <- NULL
  if (!is.null(st)) {
    sb <- .tc_smooth_bind(X, st$smooths)
    X <- sb$X; sm <- sb$smooths
  }
  list(y = as.numeric(y), X = X, formula = formula,
       terms = stats::delete.response(tt),
       xlevels = stats::.getXlevels(tt, mf),
       contrasts = ct, smooths = sm)
}

.tc_par_names <- function(X, margins, copula, W, Z, V = NULL, Q = NULL) {
  nm <- character()
  for (j in 1:3) {
    eb <- paste0("margin", j,
                 "_", .tc_family_info(margins[j])$extra_label)
    cz <- colnames(Z[[j]])
    nm <- c(nm, paste0("margin", j, "_", colnames(X[[j]])),
            if (identical(cz, "(Intercept)")) eb else paste0(eb, "_", cz))
    if (margins[j] %in% .tc_margin_names3) {
      sb <- paste0("margin", j, "_", .tc_family_info(margins[j])$shape_label)
      cq <- colnames(Q[[j]])
      nm <- c(nm, if (identical(cq, "(Intercept)")) sb else paste0(sb, "_", cq))
    }
  }
  base <- c("atanh_partial_rho12", "atanh_partial_rho13",
            "atanh_partial_rho23_given_1")
  if (copula == "pcc") {
    fam <- attr(copula, "families")
    pl <- .tc_pcc_labels(attr(copula, "root"))$long
    base <- paste0(pl, "_", fam, "_",
                   c(gaussian = "atanh_rho", student = "atanh_rho",
                     clayton = "log_theta", gumbel = "log_theta_minus_1",
                     frank = "theta", joe = "log_theta_minus_1")[.tc_pair_base(fam)])
  }
  for (k in 1:3) {
    cn <- colnames(W[[k]])
    # intercept-only predictors keep the 0.1.0 parameter names
    nm <- c(nm, if (identical(cn, "(Intercept)")) base[k] else
      paste0(base[k], "_", cn))
  }
  # intercept-only df predictors keep the earlier names
  dfn <- function(b, Vk) {
    if (.tc_df_const(Vk)) b else paste0(b, "_", colnames(Vk))
  }
  if (copula == "student") nm <- c(nm, dfn("log_df_minus_2.05", V[[1]]))
  if (copula == "pcc") {
    for (k in which(.tc_pair_base(fam) == "student")) {
      nm <- c(nm, dfn(paste0(pl[k], "_log_df_minus_2.05"), V[[k]]))
    }
  }
  nm
}

# Copula specification passed to the internal evaluators: the copula name,
# carrying the pair-copula families (and the vine root, if not 1) as
# attributes for a PCC.
.tc_cop <- function(object) {
  if (object$copula != "pcc") return(object$copula)
  root <- object$root
  cp <- if (is.null(root) || root == 1L) structure("pcc", families = object$families) else
    structure("pcc", families = object$families, root = as.integer(root))
  if (!is.null(object$nonsimplified)) attr(cp, "ns") <- object$nonsimplified
  cp
}

# Design matrices of the df predictors as used by the evaluators (NULL for
# fits without them, which gives the constant-df code path).
.tc_V <- function(object) object$df_design

# Design matrices of the three dependence predictors
# atanh(rho12), atanh(rho13), atanh(rho23|1).
.tc_make_W <- function(dependence, data, n, what = "dependence", basis = NULL) {
  if (inherits(dependence, "formula")) dependence <- rep(list(dependence), 3L)
  if (!is.list(dependence) || length(dependence) != 3L ||
      !all(vapply(dependence, inherits, logical(1), "formula"))) {
    stop(what, " must be a one-sided formula or a list of three.", call. = FALSE)
  }
  lapply(1:3, function(i) {
    f <- dependence[[i]]
    if (length(f) != 2L) stop(what, " formulas must be one-sided, e.g. ~ x.", call. = FALSE)
    st <- if (.tc_has_smooth(f)) .tc_smooth_setup(f, data, basis[[i]])
    mf <- stats::model.frame(if (is.null(st)) f else st$pf, data = data,
                             na.action = stats::na.fail)
    tt <- stats::terms(mf)
    W <- stats::model.matrix(tt, data = mf)
    if (nrow(W) != n) W <- matrix(1, n, 1, dimnames = list(NULL, "(Intercept)"))
    if (ncol(W) == 0L) stop("each ", what, " formula needs at least one term.", call. = FALSE)
    ct <- attr(W, "contrasts")
    sm <- NULL
    if (!is.null(st)) {
      sb <- .tc_smooth_bind(W, st$smooths)
      W <- sb$X; sm <- sb$smooths
    }
    list(W = W, terms = tt, xlevels = stats::.getXlevels(tt, mf),
         contrasts = ct, formula = f, smooths = sm)
  })
}

# Two-stage (inference-functions-for-margins) estimator: marginal ML for each
# equation, then copula ML with the margins held fixed. Used as the starting
# point of the joint optimization and reported for comparison.
.tc_ifm <- function(y, X, margins, copula, idx, W, Z, V = NULL, Q = NULL,
                    wt = NULL, obs = NULL) {
  # a non-simplified vine starts from the simplified one (spline coefficients 0)
  if (!is.null(attr(copula, "ns"))) attr(copula, "ns") <- NULL
  par <- numeric(idx$npar)
  U <- matrix(NA_real_, length(y[[1]]), 3)
  conv <- integer(3)
  plain <- is.null(wt) && is.null(obs)
  w1 <- if (is.null(wt)) rep(1, length(y[[1]])) else wt
  ct <- attr(y, "ct")
  for (j in 1:3) {
    wj <- if (plain) NULL else if (is.null(obs)) w1 else w1 * obs[, j]
    mj <- if (is.null(ct)) .tc_fit_margin(y[[j]], X[[j]], margins[j], Z[[j]], Q[[j]], w = wj) else
      .tc_fit_margin_ct(j, y, X, margins, Z, Q, wj, ct)
    par[idx$margin[[j]]] <- mj$par
    conv[j] <- mj$convergence
    if (!is.null(ct)) next
    if (length(idx$shapel[[j]])) {
      pj <- ncol(X[[j]]); ej <- ncol(Z[[j]])
      U[, j] <- .tc_margin_eval(y[[j]], X[[j]], margins[j], mj$par[seq_len(pj)],
                                drop(Z[[j]] %*% mj$par[pj + seq_len(ej)]),
                                shape = drop(Q[[j]] %*% mj$par[-seq_len(pj + ej)]))$u
      next
    }
    U[, j] <- .tc_margin_eval(y[[j]], X[[j]], margins[j],
                              mj$par[seq_len(ncol(X[[j]]))],
                              drop(Z[[j]] %*% mj$par[-seq_len(ncol(X[[j]]))]))$u
  }
  mc <- NULL
  if (!is.null(ct)) {
    # truncation: rescaled PITs; censoring: limits of every coordinate, and
    # the midpoints of the censored ones for the starting values
    mc <- .tc_margins_ct(par, y, X, margins, idx, Z, Q, FALSE, ct)
    U <- do.call(cbind, lapply(mc$me, `[[`, "u"))
    if (is.null(mc$cen)) mc <- NULL else U[mc$cen] <- .tc_clip((mc$Lo[mc$cen] + mc$Up[mc$cen]) / 2)
  }
  if (!is.null(obs)) U[!obs] <- NA_real_
  cop_eval <- function(th, deriv = FALSE) {
    if (is.null(mc)) return(.tc_cop_eval(th, U, W, copula, deriv = deriv, V = V))
    .tc_cop_eval_cens(th, mc$Lo, mc$Up, mc$cen, W, copula, deriv = deriv, V = V,
                      nodes = ct$nodes,
                      method = if (is.null(ct$method)) "quadrature" else ct$method)
  }
  # Start: constant dependence from the normal-score correlations placed on
  # each predictor's intercept, all other dependence coefficients zero.
  Rn <- if (is.null(obs)) stats::cor(stats::qnorm(U)) else
    stats::cor(stats::qnorm(U), use = "pairwise.complete.obs")
  c0 <- .tc_dep_inverse(Rn)
  th0 <- numeric(0)
  for (k in 1:3) {
    v <- numeric(ncol(W[[k]]))
    ic <- match("(Intercept)", colnames(W[[k]]))
    if (!is.na(ic)) v[ic] <- c0[k]
    th0 <- c(th0, v)
  }
  if (copula == "student") {
    if (.tc_df_const(V[[1]])) th0 <- c(th0, log(8 - 2.05)) else {
      v <- numeric(ncol(V[[1]]))
      ic <- match("(Intercept)", colnames(V[[1]]))
      if (!is.na(ic)) v[ic] <- log(8 - 2.05)
      th0 <- c(th0, v)
    }
  }
  if (copula == "pcc") {
    th0 <- .tc_pcc_start(U[, .tc_pcc_order(attr(copula, "root")), drop = FALSE], W,
                         attr(copula, "families"), V)
  }
  fn <- function(th) {
    v <- try(suppressWarnings({
      l <- cop_eval(th)$lc
      if (is.null(wt)) -sum(l) else -sum(wt * l)
    }), silent = TRUE)
    if (inherits(v, "try-error") || !is.finite(v)) 1e100 else v
  }
  gr <- function(th) {
    z <- try(suppressWarnings(cop_eval(th, TRUE)), silent = TRUE)
    if (inherits(z, "try-error") || is.null(z$score)) return(rep(0, length(th)))
    g <- if (is.null(wt)) -colSums(z$score) else -colSums(wt * z$score)
    g[!is.finite(g)] <- 0
    g
  }
  opt <- stats::optim(th0, fn, gr, method = "BFGS",
                      control = list(maxit = 500, reltol = 1e-10))
  par[c(idx$dep, idx$df)] <- opt$par
  list(par = par, convergence = c(conv, opt$convergence))
}

# Newton refinement of a quasi-Newton optimum, using the Jacobian of the
# analytic gradient. Steps are accepted only if they do not decrease the
# log-likelihood.
.tc_polish <- function(par, nll, ngr, maxit = 5L, tol = 1e-6) {
  val <- nll(par)
  iter <- 0L
  for (k in seq_len(maxit)) {
    g <- ngr(par)
    if (max(abs(g)) < tol) break
    H <- try(numDeriv::jacobian(ngr, par, method = "simple"), silent = TRUE)
    if (inherits(H, "try-error") || any(!is.finite(H))) break
    step <- try(solve((H + t(H)) / 2, g), silent = TRUE)
    if (inherits(step, "try-error") || any(!is.finite(step))) break
    lam <- 1
    improved <- FALSE
    for (h in 1:10) {
      cand <- par - lam * step
      cval <- nll(cand)
      if (is.finite(cval) && cval <= val) {
        par <- cand; val <- cval; improved <- TRUE
        break
      }
      lam <- lam / 2
    }
    if (!improved) break
    iter <- iter + 1L
  }
  list(par = par, value = val, iterations = iter)
}

#' Fit a trivariate continuous copula regression by full maximum likelihood
#'
#' Jointly estimates three marginal regression equations and a trivariate
#' Gaussian or Student-t copula by maximizing the full log-likelihood.
#'
#' @details
#' **Margins.** Each margin has a linear predictor `eta = X beta` and one
#' additional parameter estimated on the log scale: `gaussian`
#' (`y ~ N(eta, sigma^2)`), `lognormal` (`log y ~ N(eta, sigma^2)`) and
#' `gamma` (mean `exp(eta)`, shape `exp(extra)`), `weibull` (scale
#' `exp(eta)`, shape `exp(extra)`; the mean is `scale * gamma(1 + 1/shape)`),
#' `invgauss` (mean `mu = exp(eta)`, `sigma = exp(extra)`, variance
#' `mu^3 sigma^2`) and `beta` for responses in (0, 1) (mean `plogis(eta)`,
#' precision `phi = exp(extra)`, i.e. `Beta(mu phi, (1 - mu) phi)`),
#' `gpareto` (generalized Pareto for exceedances `y > 0` over a known
#' threshold that the user has subtracted: scale `sigma = exp(eta)`, shape
#' `xi = exp(extra) - 0.5 > -0.5`, `F(y) = 1 - (1 + xi y / sigma)^(-1 / xi)`),
#' and the three-parameter `dagum` and `singhmaddala` (scale `b = exp(eta)`,
#' `a = exp(extra)` and a third parameter with its own predictor, see
#' `shape`; all scores are analytic). For the
#' inverse Gaussian and Beta margins the derivatives of the CDF with respect
#' to the parameters use central differences. The additional parameter
#' may depend on covariates through `scale`, e.g.
#' `log sigma_i = z_i' gamma`, which gives a location-scale regression in
#' each margin. Summary functions then report the sample-average scale or
#' shape.
#'
#' **Dependence.** The correlation matrix is parameterized by
#' `atanh(rho12)`, `atanh(rho13)` and `atanh(rho23|1)`, where `rho23|1` is
#' the partial correlation, so every real parameter vector yields a
#' positive-definite matrix. The Student-t degrees of freedom are
#' `2.05 + exp(.)`. Each of the three dependence parameters may depend on
#' covariates through `dependence`, e.g. `atanh(rho12_i) = w_i' delta_12`;
#' the correlation matrix is then observation-specific and remains positive
#' definite for every observation. In that case `fit$rho` holds the fitted
#' correlations per observation, and `fit$correlation`, [tc_dependence()],
#' [tc_tau()] and [tc_taildep()] refer to their sample averages. The degrees
#' of freedom are constant unless `df` gives a predictor for
#' `log(df - 2.05)`.
#'
#' **Pair-copula construction.** With `copula = "pcc"` the copula density is
#' `c12(u1, u2) c13(u1, u3) c23|1(h12(u2 | u1), h13(u3 | u1))`, a simplified
#' C-vine, estimated jointly with the margins. With three Gaussian pairs it
#' equals the trivariate Gaussian copula in the parameterization above. The
#' derivatives of each pair copula's log-density and h-function with respect
#' to its two arguments and its parameter use central differences; they are
#' combined analytically by the chain rule. For a PCC fit, `fit$pair` and
#' `fit$ktau` hold the pair parameters and Kendall's tau per observation, and
#' `fit$correlation` is `NULL`.
#'
#' **Smooth terms.** Every formula (means, `scale`, `shape`, `dependence`,
#' `df`) may contain penalized smooths `s(x)` in the notation of
#' \pkg{mgcv}, e.g. `y ~ z + s(x)` or `dependence = ~ s(x, k = 6)`. Bases and
#' penalties come from [mgcv::smoothCon()] with the identifiability
#' constraint absorbed. `s()`, `te()` and `ti()` are supported, with or
#' without a `by` variable, including smooths with several penalties
#' (`t2()` is not). The penalized log-likelihood
#' `l - 1/2 sum_k lambda_k b_k' S_k b_k` is maximized, and the smoothing
#' parameters are chosen by the generalized Fellner-Schall iteration of
#' Wood and Fasiolo (2017) unless `sp` fixes them. For such fits `vcov()` is
#' the Bayesian covariance matrix `(H + S_lambda)^-1`, `fit$edf` holds the
#' effective degrees of freedom per coefficient, `logLik()`, `AIC()` and
#' `BIC()` use their sum, and `fit$logLik` is the unpenalized
#' log-likelihood. [tc_smooth()] returns a fitted smooth with pointwise
#' intervals. The two-stage start is unpenalized. Unpenalized regression
#' splines such as `splines::ns(x, 4)` work in every formula without any of
#' this.
#'
#' **Estimation.** Unless `start` is supplied, starting values come from a
#' two-stage (IFM) fit: marginal ML for each equation, then copula ML given
#' the margins. The joint log-likelihood is then maximized by BFGS with an
#' analytic gradient (two scalar terms, the Gamma CDF shape derivative and
#' the Student-t degrees-of-freedom derivative, use central differences),
#' optionally from several starts, and refined by Newton steps. The observed
#' information is the numerical Jacobian of the analytic gradient.
#'
#' @param formulas List of three one-response formulas.
#' @param data Data frame containing all responses and predictors. Rows with
#'   missing values are an error unless `missing` says otherwise.
#' @param weights Optional vector of non-negative prior weights, one per row
#'   of `data`. The log-likelihood is `sum(weights * l_i)`. Integer weights
#'   give the estimates, log-likelihood and model-based covariance of the
#'   data with each row repeated that many times. The sandwich covariance
#'   treats every row as one independent unit (`H^-1 sum(w_i^2 s_i s_i')
#'   H^-1`), as for sampling weights. `nobs()` and `BIC()` use the number of
#'   rows.
#' @param weights_type What the weights are. `"frequency"`: numbers of
#'   replicates; the model-based covariance is then exact and is the
#'   default of `vcov()`, `summary()`, `confint()`, [tc_dependence()] and
#'   [tc_effects()]. `"sampling"`: inverse selection probabilities or other
#'   weights that do not count observations; the model-based covariance is
#'   then too small (in the package's Monte Carlo its standard errors were
#'   51 to 84 percent of the true ones) and the sandwich covariance becomes
#'   the default of these functions. `"auto"` (the default): `"frequency"`
#'   if all weights are whole numbers, `"sampling"` otherwise. Sampling
#'   weights that happen to be whole numbers must therefore be declared.
#'   `vcov(fit, type = "hessian")` and `type = "sandwich"` always give the
#'   named matrix; `fit$vcov` is the model-based one. Ignored without
#'   `weights`.
#' @param missing What to do with missing values. `"fail"` (the default): an
#'   error. `"complete"`: complete-case analysis. `"likelihood"`: rows with
#'   a missing covariate, or with all three responses missing, are dropped;
#'   a row with one or two missing responses contributes the densities of
#'   its observed responses and the copula of those responses (the
#'   observed-data likelihood, valid if the responses are missing at
#'   random). For the Gaussian and Student-t copulas the bivariate margin is
#'   the same family; for a PCC it is a pair copula if the root is observed
#'   and otherwise a one-dimensional numerical integral over the root
#'   (Gauss-Legendre, 128 nodes), which makes such rows about 130 times more
#'   expensive. `fit$observed` is the matrix of observed responses,
#'   `fit$rows_dropped` the number of rows removed, `fit$data` the rows used.
#' @param margins Three names from `"gaussian"`, `"lognormal"`, `"gamma"`,
#'   `"weibull"`, `"invgauss"`, `"beta"`, `"gpareto"`, `"dagum"`,
#'   `"singhmaddala"`.
#' @param copula `"gaussian"`, `"student"` (trivariate elliptical copulas) or
#'   `"pcc"` (pair-copula construction, see `families`).
#' @param families For `copula = "pcc"`: three pair-copula families for
#'   `c12`, `c13` and `c23|1` of a C-vine rooted at the first response, each
#'   one of `"gaussian"`, `"student"`, `"clayton"`, `"gumbel"`, or a rotated
#'   Clayton/Gumbel (`"clayton90"`, `"clayton180"`, `"clayton270"`,
#'   `"gumbel90"`, `"gumbel180"`, `"gumbel270"`). Clayton has lower-tail
#'   dependence and Gumbel upper-tail dependence; the 180-degree rotations
#'   swap the tails and the 90/270-degree rotations give negative dependence.
#'   Each Student-t pair has its own degrees of freedom. `dependence` then
#'   gives the predictors of `atanh(rho)` (Gaussian, Student-t), `log(theta)`
#'   (Clayton) or `log(theta - 1)` (Gumbel) of the three pairs. The ordering
#'   of the responses matters: the first response is the root of the vine
#'   unless `root` is given. Further families: `"frank"` (parameter `theta`
#'   itself, any real number, negative values give negative dependence; no
#'   tail dependence) and `"joe"` with rotations `"joe90"`, `"joe180"`,
#'   `"joe270"` (`log(theta - 1)`; upper-tail dependence).
#' @param shape One-sided formula, or a list of three, for the third
#'   parameter of the three-parameter margins: `log p` of `"dagum"`
#'   (`F(y) = (1 + (y / b)^(-a))^(-p)`) and `log q` of `"singhmaddala"`
#'   (`F(y) = 1 - (1 + (y / b)^a)^(-q)`), where `b = exp(eta)` and
#'   `log a` is the `scale` predictor. Entries for two-parameter margins are
#'   ignored. The default `~ 1` gives a constant third parameter.
#' @param sp Optional vector of smoothing parameters, one per penalty in
#'   the order of `fit$sp`; if supplied they are held fixed.
#' @param smooth_basis Optional `fit$smooths` of an earlier fit of the same
#'   formulas: its spline bases (knots and constraints) are then reused
#'   instead of being rebuilt from `data`. Used by [tc_bootstrap()].
#' @param df One-sided formula for `log(df - 2.05)` of the Student-t copula,
#'   or of every Student-t pair of a PCC (a list of three formulas gives one
#'   per pair; entries of other pairs are ignored). The default `~ 1` is a
#'   constant df. With covariates, `fit$df` is the sample average and
#'   `fit$df_i` holds the values per observation.
#' @param root For `copula = "pcc"`: the response (1, 2 or 3) at the root of
#'   the C-vine. With `root = 2` the pairs are (1,2), (2,3) and (1,3 | 2), in
#'   this order in `families`, `dependence` and the coefficient names; with
#'   `root = 3` they are (1,3), (2,3) and (1,2 | 3). The root is the first
#'   argument of the two unconditional pair copulas, which matters for the
#'   90- and 270-degree rotations only. `fit$rosenblatt` is in vine order.
#' @param simplified For `copula = "pcc"`: `TRUE` (the default) fits the
#'   simplified vine, in which the conditional pair copula does not depend on
#'   the value of the conditioning variable. `FALSE` adds a penalized spline
#'   `f(u_r)` in the probability integral transform `u_r = F_r(y_r | x)` of
#'   the root to the predictor of the conditional pair
#'   (`eta_i = w_i' delta + f(u_r,i)`, with `f` integrating to zero over
#'   (0, 1)). All parameters are still estimated jointly; `u_r` moves with
#'   the parameters of the root's margin. The smoothing parameter is
#'   estimated like those of `s()` terms (it is the last element of `sp`),
#'   the smooth is listed in `summary()` and returned by [tc_smooth()] (as a
#'   function of `u`), and [tc_simplifying_test()] tests the simplified vine
#'   against it. [tc_dependence()], [tc_tau()] and `fit$pair` then refer to
#'   the fitted conditioning values of the sample.
#' @param k_cond Number of cubic B-spline basis functions of `f` (at least
#'   4; `f` has `k_cond - 1` coefficients). Ignored if `simplified = TRUE`.
#'   The default is 10 since 0.17.0 (6 before; `k_cond = 6` gives the fits
#'   of 0.15.0 and 0.16.0): with 6 the pointwise intervals of a curve with a
#'   turning point covered too rarely.
#' @param truncation Optional list of three, one entry per response: `NULL`
#'   (not truncated), `c(lower, upper)`, or `list(lower = , upper = )` where
#'   each bound is a number or the name of a column of `data` (bounds that
#'   differ between observations, such as entry times). Use `-Inf` or `Inf`
#'   for a side without bound. All values of the response must lie inside.
#' @param truncation_type `"margin"` (the default): each truncated response
#'   has the truncated distribution as its margin, with density
#'   `f(y) / (F(b) - F(a))` and probability integral transform
#'   `(F(y) - F(a)) / (F(b) - F(a))`, and the copula joins these margins; the
#'   model then describes the truncated population. `"joint"`: the three
#'   responses follow the untruncated model and a row is observed only if
#'   the truncated response lies inside its bounds; the likelihood is the
#'   untruncated joint density divided by `F(b) - F(a)`, and the parameters
#'   describe the untruncated population. `"joint"` is implemented for one
#'   truncated response.
#' @param censoring Optional list of three, one entry per response: `NULL`
#'   (never censored) or `c(status = "column", upper = "column")`. The
#'   status column is coded as in `survival::Surv()`: 1 exact, 0
#'   right-censored at the response value, 2 left-censored at the response
#'   value, 3 interval-censored between the response value and the value in
#'   the `upper` column (needed only if code 3 occurs). A row contributes
#'   the densities of its exact responses and the copula density integrated
#'   over the censored ones (Gauss-Legendre quadrature; one, two or three
#'   censored responses per row, every copula). `fit$uniforms` and
#'   `fit$rosenblatt` are `NA` where a response is censored. Rows with three
#'   censored responses are slow (4096 evaluations of the copula density
#'   each). Not available together with missing responses.
#' @param cens_nodes Number of quadrature nodes per censored response for
#'   rows with one, two and three censored responses. The default
#'   `c(64, 32, 16)` gives errors of the log-likelihood contribution below
#'   1e-6 for moderate dependence (Kendall's tau up to about 0.5); with
#'   very strong dependence more nodes are needed (see `cens_check`). With
#'   `cens_method = "closed"` the three numbers are the nodes per response
#'   when one, two or three responses are left for quadrature.
#' @param cens_method How the copula density is integrated over censored
#'   responses. `"closed"` (the default since 0.17.1): one censored
#'   response per row is integrated in closed form and quadrature runs over
#'   the remaining ones, so a row with one censored response needs no
#'   quadrature, a row with two needs `cens_nodes[1]` nodes and a row with
#'   three `cens_nodes[2]^2`. `"quadrature"` (the method of 0.16.0, and the
#'   default of 0.17.0): Gauss-Legendre quadrature over every censored
#'   response; slower and less accurate, kept for reproduction. Closed forms
#'   exist for the Gaussian and Student-t copulas and, in a vine, for every
#'   response except the root; a row in which only the root is censored is
#'   computed by quadrature. The two methods give the same likelihood up to
#'   the quadrature error.
#' @param cens_check Check the accuracy of the quadrature after the fit? See
#'   [tc_cens_check()]. The result is stored as `fit$cens_accuracy`; a
#'   warning is given if an estimate would move by more than 0.1 standard
#'   errors with twice the nodes. The estimates are not affected.
#' @param start Optional starting parameter vector (unconstrained scale). If
#'   supplied, the IFM stage is skipped.
#' @param control List passed to [stats::optim()].
#' @param hessian Compute the observed information and covariance matrix?
#' @param nstart Number of optimizer starts. Starts beyond the first perturb
#'   the dependence parameters of the first start with `N(0, 0.5^2)` noise
#'   (uses the R random number generator).
#' @param polish Refine the BFGS solution with Newton steps?
#' @param dependence One-sided formula, or a list of three, for the linear
#'   predictors of `atanh(rho12)`, `atanh(rho13)` and `atanh(rho23|1)`.
#'   The default `~ 1` gives constant dependence. A single formula is used
#'   for all three parameters.
#' @param scale One-sided formula, or a list of three, for the linear
#'   predictor of the log scale (`log sigma` for Gaussian and log-normal
#'   margins) or log shape (Gamma margin) of each equation. The default
#'   `~ 1` gives constant scale and shape. A single formula is used for all
#'   three margins.
#' @return An object of class `tricopula`: a list with, among others,
#'   `coefficients` (unconstrained scale), `logLik`, `vcov`, `hessian`
#'   (of the negative log-likelihood), `gradient`, `scores`, `correlation`,
#'   `df`, `uniforms` (PITs), `ifm` (two-stage estimates) and `starts`.
#' @examples
#' set.seed(1)
#' n <- 300
#' d <- data.frame(x = rnorm(n))
#' U <- tc_rcopula(n, matrix(c(1, .5, .3, .5, 1, .4, .3, .4, 1), 3))
#' d$y1 <- qnorm(U[, 1], 1 + 0.5 * d$x, 1)
#' d$y2 <- qlnorm(U[, 2], 0.2 + 0.3 * d$x, 0.5)
#' d$y3 <- qgamma(U[, 3], shape = 2.5, scale = exp(0.1 + 0.2 * d$x) / 2.5)
#' fit <- tricopula(list(y1 ~ x, y2 ~ x, y3 ~ x), d,
#'                  margins = c("gaussian", "lognormal", "gamma"))
#' summary(fit)
#' @export
tricopula <- function(formulas, data, margins = c("gaussian", "lognormal", "gamma"),
                      copula = c("gaussian", "student", "pcc"), start = NULL,
                      control = list(maxit = 2000, reltol = 1e-9),
                      hessian = TRUE, nstart = 1L, polish = TRUE,
                      dependence = ~ 1, scale = ~ 1, families = NULL,
                      df = ~ 1, root = 1, shape = ~ 1, sp = NULL, smooth_basis = NULL,
                      weights = NULL, missing = c("fail", "likelihood", "complete"),
                      weights_type = c("auto", "frequency", "sampling"),
                      simplified = TRUE, k_cond = 10L,
                      truncation = NULL, truncation_type = c("margin", "joint"),
                      censoring = NULL, cens_nodes = c(64L, 32L, 16L),
                      cens_method = c("closed", "quadrature"), cens_check = TRUE) {
  copula <- match.arg(copula)
  cens_method <- match.arg(cens_method)
  truncation_type <- match.arg(truncation_type)
  if (!is.logical(simplified) || length(simplified) != 1L || is.na(simplified)) {
    stop("simplified must be TRUE or FALSE.")
  }
  if (!simplified && copula != "pcc") {
    stop("simplified = FALSE is only used with copula = \"pcc\".")
  }
  ns <- if (simplified) NULL else .tc_ns_setup(k_cond)
  weights_type <- match.arg(weights_type)
  if (copula == "pcc") {
    if (length(families) != 3L || any(!families %in% .tc_pair_families)) {
      stop("For copula = \"pcc\" supply families: three of ",
           paste(.tc_pair_families, collapse = ", "), ".")
    }
  } else if (!is.null(families)) {
    stop("families is only used with copula = \"pcc\".")
  }
  if (!is.numeric(root) || length(root) != 1L || !root %in% 1:3) {
    stop("root must be 1, 2 or 3.")
  }
  root <- as.integer(root)
  if (root != 1L && copula != "pcc") stop("root is only used with copula = \"pcc\".")
  cop <- if (copula != "pcc") copula else if (root == 1L) {
    structure(copula, families = families)
  } else structure(copula, families = families, root = root)
  if (!is.null(ns)) attr(cop, "ns") <- ns
  if (length(formulas) != 3L) stop("Supply exactly three formulas.")
  if (length(margins) != 3L || any(!margins %in% .tc_margin_names)) {
    stop("margins must contain exactly three supported margin names.")
  }
  missing <- match.arg(missing)
  # Prior weights and missing values (0.13.0). With weights = NULL and
  # missing = "fail" nothing below changes the data, and wt and obs stay
  # NULL: the code path of the earlier versions.
  wt <- NULL; obs <- NULL; rows_dropped <- 0L
  if (!is.null(weights)) {
    if (!is.numeric(weights) || length(weights) != nrow(data) || anyNA(weights) ||
        any(!is.finite(weights)) || any(weights < 0)) {
      stop("weights must be a vector of non-negative numbers, one per row of data.")
    }
    wt <- as.numeric(weights)
    if (weights_type == "auto") {
      weights_type <- if (all(abs(wt - round(wt)) < 1e-8)) "frequency" else "sampling"
    }
  } else weights_type <- NULL
  if (missing != "fail") {
    fl <- function(f) if (inherits(f, "formula")) list(f) else f
    vars_y <- lapply(formulas, function(f) intersect(all.vars(f[[2]]), names(data)))
    vars_all <- unique(unlist(lapply(c(formulas, fl(dependence), fl(scale), fl(df), fl(shape)),
                                     all.vars)))
    vars_x <- setdiff(intersect(vars_all, names(data)), unlist(vars_y))
    keep <- .tc_missing_rows(data, vars_y, vars_x, missing)
    rows_dropped <- sum(!keep)
    if (rows_dropped) {
      data <- data[keep, , drop = FALSE]
      if (!is.null(wt)) wt <- wt[keep]
    }
    if (!nrow(data)) stop("No rows are left after removing missing values.")
  }
  mm <- lapply(1:3, function(j) {
    .tc_make_X(formulas[[j]], data, smooth_basis$mean[[j]], na_ok = missing == "likelihood")
  })
  n <- vapply(mm, function(x) length(x$y), integer(1))
  if (length(unique(n)) != 1L) stop("All three margins must have the same number of observations.")
  y <- lapply(mm, `[[`, "y")
  X <- lapply(mm, `[[`, "X")
  if (missing == "likelihood" && any(vapply(y, anyNA, logical(1)))) {
    obs <- !vapply(y, is.na, logical(n[1]))
    obs <- matrix(obs, n[1], 3)
    # placeholder values keep the evaluators vectorized; .tc_eval removes
    # every contribution of a missing response
    for (j in 1:3) y[[j]][!obs[, j]] <- stats::median(y[[j]][obs[, j]])
  }
  for (j in 1:3) .tc_check_support(y[[j]], margins[j])
  # Truncation and censoring (0.16.0). With both arguments NULL `ct` is NULL
  # and y carries no attribute: the code path of the earlier versions.
  ct <- NULL
  if (!is.null(truncation) || !is.null(censoring)) {
    if (!is.null(censoring) && !is.null(obs)) {
      stop("censoring cannot be combined with missing responses (missing = \"likelihood\").")
    }
    if (!is.numeric(cens_nodes) || length(cens_nodes) != 3L || anyNA(cens_nodes) ||
        any(cens_nodes < 4) || any(cens_nodes > 400)) {
      stop("cens_nodes must be three numbers between 4 and 400.")
    }
    ct <- .tc_ct_setup(truncation, truncation_type, censoring, data, y, margins)
    if (!is.null(ct)) ct$nodes <- as.integer(cens_nodes)
    if (!is.null(ct)) ct$method <- cens_method
    if (!is.null(ct)) attr(y, "ct") <- ct
  }
  p <- vapply(X, ncol, integer(1))
  dd <- .tc_make_W(dependence, data, n[1], basis = smooth_basis$dep)
  W <- lapply(dd, `[[`, "W")
  zz <- .tc_make_W(scale, data, n[1], what = "scale", basis = smooth_basis$scale)
  Z <- lapply(zz, `[[`, "W")
  # third (shape) predictor of the three-parameter margins
  qq <- .tc_make_W(shape, data, n[1], what = "shape", basis = smooth_basis$shape)
  is3 <- margins %in% .tc_margin_names3
  if (!any(is3) && any(!vapply(qq, function(v) .tc_df_const(v$W), logical(1)))) {
    stop("shape is only used with three-parameter margins (",
         paste(.tc_margin_names3, collapse = ", "), ").")
  }
  Q <- if (any(is3)) lapply(qq, `[[`, "W") else NULL
  # df predictors: one for the Student-t copula (the first formula), one per
  # Student-t pair of a PCC. V stays NULL when all of them are intercept-only,
  # which is the constant-df code path of the earlier versions.
  ff <- .tc_make_W(df, data, n[1], what = "df", basis = smooth_basis$df)
  has_df <- if (copula == "student") c(TRUE, FALSE, FALSE) else if (copula == "pcc") {
    .tc_pair_base(families) == "student"
  } else rep(FALSE, 3L)
  df_const <- vapply(ff, function(v) .tc_df_const(v$W), logical(1))
  if (any(!df_const)) {
    if (!any(has_df)) {
      stop("df is only used with a Student-t copula or Student-t pair copulas.")
    }
    if (copula == "student" && !inherits(df, "formula")) {
      stop("For copula = \"student\" df must be a single one-sided formula.")
    }
  }
  Vd <- if (all(df_const[has_df])) NULL else if (copula == "student") {
    list(ff[[1]]$W)
  } else lapply(ff, `[[`, "W")
  idx <- .tc_index(p, cop, vapply(W, ncol, integer(1)),
                   vapply(Z, ncol, integer(1)),
                   if (is.null(Vd)) c(1L, 1L, 1L) else
                     vapply(ff, function(v) ncol(v$W), integer(1)),
                   if (is.null(Q)) c(0L, 0L, 0L) else
                     ifelse(is3, vapply(Q, ncol, integer(1)), 0L),
                   if (is.null(ns)) 0L else ns$ncoef)
  par_names <- .tc_par_names(X, margins, cop, W, Z, Vd, Q)
  if (!is.null(ns)) {
    par_names <- c(par_names, paste0(.tc_pcc_labels(root)$long[3], "_s(u", root, ").",
                                     seq_len(ns$ncoef)))
  }

  nll <- function(par) {
    # Line-search excursions to extreme parameters may produce NaN; they are
    # rejected through the penalty value, so the warnings carry no information.
    z <- try(suppressWarnings(.tc_eval(par, y, X, margins, cop, idx,
                                       W = W, Z = Z, V = Vd, Q = Q, wt = wt,
                                       obs = obs)$loglik),
             silent = TRUE)
    if (inherits(z, "try-error") || !is.finite(z)) return(1e100)
    -z
  }
  ngr <- function(par) {
    z <- try(suppressWarnings(.tc_eval(par, y, X, margins, cop, idx,
                                       deriv = TRUE, W = W, Z = Z, V = Vd, Q = Q,
                                       wt = wt, obs = obs)),
             silent = TRUE)
    if (inherits(z, "try-error") || is.null(z$score)) return(rep(0, idx$npar))
    g <- -colSums(z$score)
    g[!is.finite(g)] <- 0
    g
  }

  ifm <- NULL
  if (is.null(start)) {
    ifm <- .tc_ifm(y, X, margins, cop, idx, W, Z, Vd, Q, wt, obs)
    init <- ifm$par
    ifm <- list(coefficients = stats::setNames(ifm$par, par_names),
                logLik = -nll(ifm$par), convergence = ifm$convergence)
  } else {
    if (length(start) != idx$npar) stop("start has incorrect length.")
    init <- as.numeric(start)
  }

  # Penalized smooth terms (0.12.0). Without s() terms `pens` is empty and
  # everything below is the unpenalized code of the earlier versions.
  pl_dep <- if (copula == "pcc") .tc_pcc_labels(root)$short else c("12", "13", "23|1")
  blocks <- list()
  for (j in 1:3) {
    fi <- .tc_family_info(margins[j])
    blocks <- c(blocks, list(
      list(prefix = paste0("margin", j), positions = idx$beta[[j]],
           smooths = mm[[j]]$smooths),
      list(prefix = paste0("margin", j, "_", fi$extra_label),
           positions = idx$extral[[j]], smooths = zz[[j]]$smooths)))
    if (is3[j]) {
      blocks <- c(blocks, list(list(prefix = paste0("margin", j, "_", fi$shape_label),
                                    positions = idx$shapel[[j]],
                                    smooths = qq[[j]]$smooths)))
    }
  }
  for (k in 1:3) {
    blocks <- c(blocks, list(list(prefix = paste0("dep", pl_dep[k]),
                                  positions = idx$depl[[k]], smooths = dd[[k]]$smooths)))
  }
  if (!is.null(Vd)) {
    for (k in seq_along(idx$dfl)) {
      if (!length(idx$dfl[[k]])) next
      blocks <- c(blocks, list(list(
        prefix = if (copula == "student") "df" else paste0("df", pl_dep[k]),
        positions = idx$dfl[[k]], smooths = ff[[k]]$smooths)))
    }
  }
  pc <- .tc_collect_penalties(blocks)
  pens <- pc$pens
  sblocks <- pc$sblocks
  if (!is.null(ns)) {
    # the spline in the conditioning value: the last smooth and penalty
    lab <- paste0("dep", pl_dep[3], ":s(u", root, ")")
    sblocks <- c(sblocks, list(list(label = lab, ind = idx$ns,
                                    smooth = list(tc_ns = TRUE, label = paste0("s(u", root, ")")))))
    pens <- c(pens, list(list(ind = idx$ns, S = ns$S, block = length(sblocks), label = lab)))
  }
  lam <- NULL
  sp_iterations <- 0L; sp_converged <- NA
  if (length(pens)) {
    if (is.null(sp)) lam <- rep(1, length(pens)) else {
      if (!is.numeric(sp) || length(sp) != length(pens) || any(sp < 0)) {
        stop("sp must be a vector of ", length(pens), " non-negative smoothing parameters (",
             paste(vapply(pens, `[[`, "", "label"), collapse = ", "), ").")
      }
      lam <- as.numeric(sp)
    }
    nll_u <- nll; ngr_u <- ngr
    nll <- function(par) nll_u(par) + .tc_pen_value(par, pens, lam)
    ngr <- function(par) ngr_u(par) + .tc_pen_grad(par, pens, lam)
    # The start is unpenalized, and a large penalty makes the problem stiff
    # for BFGS: take Newton steps on the penalized objective first.
    init <- .tc_polish(init, nll, ngr, maxit = 10L)$par
  } else if (!is.null(sp)) {
    stop("sp is only used with s() terms or simplified = FALSE.")
  }

  nstart <- max(1L, as.integer(nstart))
  inits <- vector("list", nstart)
  inits[[1]] <- init
  for (s in seq_len(nstart - 1L)) {
    cand <- init
    for (k in 1:3) {
      ii <- idx$depl[[k]]
      cand[ii] <- cand[ii] + stats::rnorm(length(ii), sd = c(0.5, rep(0.1, length(ii) - 1L)))
    }
    inits[[s + 1L]] <- cand
  }
  opts <- lapply(inits, function(st) {
    stats::optim(st, nll, ngr, method = "BFGS", control = control)
  })
  values <- vapply(opts, `[[`, numeric(1), "value")
  opt <- opts[[which.min(values)]]
  par <- opt$par
  newton <- 0L
  # (when smoothing parameters are still to be selected, the refinement is
  # left to the selection loop below)
  if (polish && !(length(pens) && is.null(sp))) {
    pol <- .tc_polish(par, nll, ngr)
    par <- pol$par
    newton <- pol$iterations
  }
  if (length(pens) && is.null(sp)) {
    # smoothing-parameter selection: generalized Fellner-Schall iteration,
    # each inner fit started from the previous estimate
    # Stops when every smoothing parameter changes by less than 1% or the
    # effective degrees of freedom of every smooth by less than 0.002.
    # The plain update is slow for a smoothing parameter on its way to
    # infinity (an effect that is really in the penalty's null space): the
    # log-step of a parameter that keeps moving in the same direction is
    # therefore doubled at each iteration (up to 64 times the plain step and
    # a factor of 1000 per iteration) and reset when the direction changes.
    # The numerical information matrix is the expensive part: after the
    # first three iterations it is recomputed only when a coefficient has
    # moved by more than 5% (relative to 0.1 + its size) since it was last
    # computed, and convergence is only declared with a fresh one.
    sp_converged <- FALSE
    edf_old <- NULL
    Hu <- NULL; par_H <- NULL
    acc <- rep(1, length(pens)); dir_old <- rep(0, length(pens))
    for (it in 1:80) {
      fresh <- is.null(Hu) || it <= 3L ||
        max(abs(par - par_H) / (0.1 + abs(par_H))) > 0.05
      if (fresh) {
        Hu <- try(numDeriv::jacobian(ngr_u, par, method = "simple"), silent = TRUE)
        if (inherits(Hu, "try-error") || any(!is.finite(Hu))) break
        Hu <- (Hu + t(Hu)) / 2
        par_H <- par
      }
      up <- .tc_fs_update(par, Hu, pens, lam, length(sblocks))
      if (is.null(up)) break
      sp_iterations <- it
      step <- log(up$lam) - log(lam)
      if (max(abs(step)) < 0.01 ||
          (!is.null(edf_old) && max(abs(up$edf - edf_old)) < 0.002)) {
        if (fresh) {
          sp_converged <- TRUE
          break
        }
        Hu <- NULL                         # confirm with a fresh matrix
        next
      }
      edf_old <- up$edf
      same <- sign(step) == dir_old & abs(step) > 0.01
      acc <- ifelse(same, pmin(2 * acc, 64), 1)
      dir_old <- sign(step)
      lam <- pmin(pmax(lam * exp(pmin(pmax(step * acc, -log(1e3)), log(1e3))), 1e-7), 1e7)
      # Inner fit from the previous estimate. First Newton steps that reuse
      # the information matrix just computed (one gradient per step, no new
      # numerical Hessian); full Newton steps and then BFGS only if those
      # do not bring the gradient down.
      nf <- .tc_newton_fixed(par, nll, ngr, Hu + .tc_pen_matrix(idx$npar, pens, lam))
      par <- nf$par
      if (!nf$ok) {
        pol <- .tc_polish(par, nll, ngr, maxit = 10L)
        par <- pol$par
        newton <- pol$iterations
      }
      if (max(abs(ngr(par))) > 1e-4) {
        opt <- stats::optim(par, nll, ngr, method = "BFGS", control = control)
        par <- opt$par
        if (polish) {
          pol <- .tc_polish(par, nll, ngr)
          par <- pol$par
          newton <- pol$iterations
        }
      }
    }
    if (!sp_converged) {
      warning("The smoothing-parameter iteration did not converge in ", sp_iterations,
              " steps.", call. = FALSE)
    }
  }
  names(par) <- par_names

  ev <- .tc_eval(par, y, X, margins, cop, idx, deriv = TRUE, W = W, Z = Z, V = Vd, Q = Q,
                 wt = wt, obs = obs)
  if (!is.finite(ev$loglik)) stop("The log-likelihood is not finite at the optimum.")
  H <- V <- NULL
  df_boundary <- FALSE
  if (hessian) {
    # penalized fits: two Richardson iterations instead of four (half the
    # gradient evaluations; the error stays far below that of the estimates)
    H <- if (length(pens)) {
      numDeriv::jacobian(ngr, par, method.args = list(r = 2))
    } else numDeriv::jacobian(ngr, par)
    H <- (H + t(H)) / 2
    dimnames(H) <- list(par_names, par_names)
    ok <- all(is.finite(H)) &&
      min(eigen(H, symmetric = TRUE, only.values = TRUE)$values) > 0
    if (ok) {
      V <- solve(H)
      dimnames(V) <- list(par_names, par_names)
    } else if (length(idx$df) > 0L && all(is.finite(H)) &&
               min(eigen(H[-idx$df, -idx$df], symmetric = TRUE,
                         only.values = TRUE)$values) > 0) {
      # df has drifted towards infinity (Gaussian limit): the likelihood is
      # flat in the df direction. Report the covariance of the other
      # parameters conditional on df and NA for df itself.
      V <- matrix(NA_real_, idx$npar, idx$npar, dimnames = list(par_names, par_names))
      V[-idx$df, -idx$df] <- solve(H[-idx$df, -idx$df])
      df_boundary <- TRUE
    } else {
      warning("The observed information is not positive definite; ",
              "no covariance matrix is returned.", call. = FALSE)
    }
  }
  if (length(idx$df) > 0L && (df_boundary || any(ev$df > 100, na.rm = TRUE))) {
    warning("Estimated Student-t degrees of freedom are ",
            paste(signif(ev$df[!is.na(ev$df)], 3), collapse = ", "),
            ": the fit is practically a Gaussian copula and the df is not ",
            "identified", if (df_boundary) " (its standard error is NA)" else "",
            ". Consider a Gaussian copula.", call. = FALSE)
  }
  colnames(ev$score) <- par_names
  # effective degrees of freedom of a penalized fit: diag(V_p H) with H the
  # unpenalized information (the stored hessian and vcov are penalized:
  # vcov is the Bayesian covariance matrix)
  edf <- NULL; grad <- colSums(ev$score)
  if (length(pens)) {
    grad <- grad - .tc_pen_grad(par, pens, lam)
    Sl <- .tc_pen_matrix(idx$npar, pens, lam)
    Hp <- if (!is.null(H)) H else {
      Hx <- numDeriv::jacobian(ngr, par, method = "simple")
      (Hx + t(Hx)) / 2
    }
    Vp <- if (!is.null(V) && !anyNA(V)) V else try(solve(Hp), silent = TRUE)
    if (!inherits(Vp, "try-error")) {
      edf <- stats::setNames(rowSums(Vp * (Hp - Sl)), par_names)
    }
    names(lam) <- vapply(pens, `[[`, "", "label")
  }
  resp <- vapply(formulas, function(f) deparse(f[[2]]), character(1))
  fit <- structure(list(
    coefficients = par, logLik = ev$loglik, n = n[1],
    convergence = opt$convergence, message = opt$message,
    hessian = H, vcov = V, copula = copula, margins = margins,
    formulas = formulas, design = X, response = y,
    correlation = ev$R, df = ev$df, fitted_margins = ev$margins,
    families = families, pair = ev$pair, ktau = ev$ktau,
    rosenblatt = ev$rosenblatt,
    uniforms = ev$U, loglik_i = ev$ll_i, optim = opt,
    parameter_count = length(par),
    gradient = stats::setNames(grad, par_names),
    scores = ev$score, index = idx, ifm = ifm,
    starts = list(par = inits, logLik = -values), newton_iterations = newton,
    terms = lapply(mm, `[[`, "terms"), xlevels = lapply(mm, `[[`, "xlevels"),
    contrasts = lapply(mm, `[[`, "contrasts"), response_names = resp,
    data = data, control = control, df_boundary = df_boundary,
    scale_design = Z, scale_formulas = lapply(zz, `[[`, "formula"),
    scale_terms = lapply(zz, `[[`, "terms"),
    scale_xlevels = lapply(zz, `[[`, "xlevels"),
    scale_contrasts = lapply(zz, `[[`, "contrasts"),
    scale_constant = all(vapply(Z, function(w) identical(colnames(w), "(Intercept)"), logical(1))),
    rho = ev$rho, dep_design = W, dependence = lapply(dd, `[[`, "formula"),
    dep_terms = lapply(dd, `[[`, "terms"),
    dep_xlevels = lapply(dd, `[[`, "xlevels"),
    dep_contrasts = lapply(dd, `[[`, "contrasts"),
    dep_constant = all(vapply(W, function(w) identical(colnames(w), "(Intercept)"), logical(1))),
    root = root, df_i = ev$df_i, df_design = Vd, df_constant = is.null(Vd),
    df_formulas = lapply(ff, `[[`, "formula"), df_terms = lapply(ff, `[[`, "terms"),
    df_xlevels = lapply(ff, `[[`, "xlevels"),
    df_contrasts = lapply(ff, `[[`, "contrasts"),
    shape_design = Q, shape_formulas = lapply(qq, `[[`, "formula"),
    shape_terms = lapply(qq, `[[`, "terms"),
    shape_xlevels = lapply(qq, `[[`, "xlevels"),
    shape_contrasts = lapply(qq, `[[`, "contrasts"),
    penalties = if (length(pens)) pens, sp = lam, edf = edf,
    smooth_blocks = if (length(sblocks)) sblocks,
    weights = wt, weights_type = weights_type, observed = obs, missing = missing, rows_dropped = rows_dropped,
    sp_iterations = sp_iterations, sp_converged = sp_converged,
    penalized_logLik = if (length(pens)) ev$loglik - .tc_pen_value(par, pens, unname(lam)),
    smooths = list(mean = lapply(mm, `[[`, "smooths"), scale = lapply(zz, `[[`, "smooths"),
                   shape = lapply(qq, `[[`, "smooths"), dep = lapply(dd, `[[`, "smooths"),
                   df = lapply(ff, `[[`, "smooths")),
    simplified = simplified, nonsimplified = ns,
    truncation = if (!is.null(ct) && any(!vapply(ct$lower, is.null, logical(1)))) truncation,
    truncation_type = if (!is.null(ct) && any(!vapply(ct$lower, is.null, logical(1)))) truncation_type,
    censoring = if (!is.null(ct$cens)) censoring, censored = ct$cens,
    cens_nodes = if (!is.null(ct$cens)) ct$nodes,
    cens_method = if (!is.null(ct$cens)) cens_method,
    call = match.call()
  ), class = "tricopula")
  # Accuracy of the quadrature (0.17.0): one more evaluation with twice the
  # nodes; the estimates are not changed.
  if (!is.null(ct$cens) && isTRUE(cens_check)) {
    ck <- try(tc_cens_check(fit), silent = TRUE)
    if (!inherits(ck, "try-error")) {
      fit$cens_accuracy <- ck
      if (isTRUE(ck$max_shift > 0.1) || (is.na(ck$max_shift) && abs(ck$loglik_change) > 0.01)) {
        warning("The quadrature over the censored responses is not accurate enough: with ",
                "twice the nodes ",
                if (is.na(ck$max_shift)) paste0("the log-likelihood changes by ",
                                                signif(ck$loglik_change, 2)) else
                  paste0("an estimate would move by about ", signif(ck$max_shift, 2),
                         " standard errors"),
                ". Refit with larger cens_nodes",
                if (cens_method == "quadrature") " or cens_method = \"closed\"", ".",
                call. = FALSE)
      }
    }
  }
  fit
}

.tc_check <- function(object) {
  if (!inherits(object, "tricopula")) stop("object must be a tricopula fit.", call. = FALSE)
  invisible(object)
}

#' Evaluate the joint log-likelihood or its gradient at a parameter vector
#'
#' @param object A `tricopula` fit.
#' @param par Parameter vector on the unconstrained scale.
#' @return `tc_loglik()` returns the log-likelihood; `tc_gradient()` its
#'   gradient with respect to `par`.
#' @export
tc_loglik <- function(object, par = object$coefficients) {
  .tc_check(object)
  .tc_eval(par, object$response, object$design, object$margins,
           .tc_cop(object), object$index, W = object$dep_design,
           Z = object$scale_design, V = .tc_V(object), Q = object$shape_design,
           wt = object$weights, obs = object$observed)$loglik
}

#' @rdname tc_loglik
#' @export
tc_gradient <- function(object, par = object$coefficients) {
  .tc_check(object)
  ev <- .tc_eval(par, object$response, object$design, object$margins,
                 .tc_cop(object), object$index, deriv = TRUE,
                 W = object$dep_design, Z = object$scale_design, V = .tc_V(object),
                 Q = object$shape_design, wt = object$weights, obs = object$observed)
  stats::setNames(colSums(ev$score), names(object$coefficients))
}

#' Methods for tricopula fits
#'
#' @param object,x A `tricopula` fit (or its summary).
#' @param type For `vcov()`: `"hessian"` (inverse observed information; for
#'   fits with penalized smooths the Bayesian covariance matrix
#'   `(H + S_lambda)^-1`), `"sandwich"` (robust to misspecification,
#'   `H^-1 S'S H^-1`) or, for fits with penalized smooths, `"corrected"`:
#'   the Bayesian covariance matrix plus a first-order correction for the
#'   uncertainty of the smoothing parameters (Wood, Pya and Saefken, 2016).
#'   The correction refits the model once per smoothing parameter. The
#'   default is `"hessian"`, except for fits with sampling weights (see
#'   `weights_type` of [tricopula()]), where it is `"sandwich"`.
#' @param level Confidence level.
#' @param parm Optional names or indices of the parameters to report.
#' @param ... Unused.
#' @name tricopula-methods
NULL

#' @rdname tricopula-methods
#' @export
logLik.tricopula <- function(object, ...) {
  # penalized fits: effective degrees of freedom
  structure(object$logLik,
            df = if (is.null(object$edf)) object$parameter_count else sum(object$edf),
            nobs = object$n, class = "logLik")
}

#' @rdname tricopula-methods
#' @export
nobs.tricopula <- function(object, ...) object$n

#' @rdname tricopula-methods
#' @export
vcov.tricopula <- function(object, type = NULL, ...) {
  type <- .tc_vcov_type(object, type, c("hessian", "sandwich", "corrected"))
  if (type == "corrected") return(.tc_sp_correction(object))
  if (is.null(object$vcov)) stop("Hessian could not be inverted reliably.")
  if (type == "hessian") return(object$vcov)
  keep <- !is.na(diag(object$vcov))   # df excluded when at the Gaussian limit
  V <- object$vcov
  V[keep, keep] <- object$vcov[keep, keep] %*%
    crossprod(object$scores[, keep, drop = FALSE]) %*% object$vcov[keep, keep]
  V
}

# Covariance type used when the caller gives none: the sandwich for fits with
# sampling weights (0.13.1), the inverse information otherwise.
.tc_vcov_type <- function(object, type, choices = c("hessian", "sandwich")) {
  if (is.null(type)) {
    return(if (identical(object$weights_type, "sampling")) "sandwich" else "hessian")
  }
  match.arg(type, choices)
}

# Natural-scale parameters as a function of the unconstrained vector.
.tc_natural <- function(par, object) {
  idx <- object$index
  W <- object$dep_design
  sc <- vapply(1:3, function(j) {
    mean(exp(drop(object$scale_design[[j]] %*% par[idx$extral[[j]]])))
  }, numeric(1))
  mnm <- paste0("margin", 1:3, "_",
                vapply(object$margins, function(m) .tc_family_info(m)$natural,
                       character(1)))
  # third parameter (sample average) of the three-parameter margins
  j3 <- which(lengths(idx$shapel) > 0L)
  sc3 <- vapply(j3, function(j) {
    mean(exp(drop(object$shape_design[[j]] %*% par[idx$shapel[[j]]])))
  }, numeric(1))
  nm3 <- vapply(j3, function(j) {
    paste0("margin", j, "_", .tc_family_info(object$margins[j])$natural3)
  }, character(1))
  if (object$copula == "pcc") {
    # sample averages of each pair's natural parameter and Kendall's tau
    fam <- object$families
    eta <- lapply(1:3, function(k) drop(W[[k]] %*% par[idx$depl[[k]]]))
    if (!is.null(object$nonsimplified)) {
      # non-simplified vine: averages over the fitted conditioning values
      u_r <- object$uniforms[, .tc_pcc_order(object$root)[1]]
      eta[[3]] <- (eta[[3]] + drop(object$nonsimplified$basis(u_r, FALSE)$B %*%
                                     par[idx$ns]))[!is.na(u_r)]
    }
    nat <- vapply(1:3, function(k) mean(.tc_pair_natural(eta[[k]], fam[k])), numeric(1))
    kt <- vapply(1:3, function(k) mean(.tc_pair_ktau(eta[[k]], fam[k])), numeric(1))
    it <- .tc_pair_base(fam) == "student"
    pl <- .tc_pcc_labels(object$root)$short
    return(stats::setNames(
      c(sc, sc3, nat, kt, .tc_mean_df(par, object)),
      c(mnm, nm3, paste0("par", pl), paste0("tau", pl),
        if (any(it)) paste0("df", pl[it]))))
  }
  dep <- .tc_dep_vec(drop(W[[1]] %*% par[idx$depl[[1]]]),
                     drop(W[[2]] %*% par[idx$depl[[2]]]),
                     drop(W[[3]] %*% par[idx$depl[[3]]]))
  out <- c(sc, sc3, mean(dep$r12), mean(dep$r13), mean(dep$r23),
           mean(dep$p23))
  nm <- c(mnm, nm3, "rho12", "rho13", "rho23", "rho23|1")
  if (object$copula == "student") {
    out <- c(out, .tc_mean_df(par, object))
    nm <- c(nm, "df")
  }
  stats::setNames(out, nm)
}

# Degrees of freedom implied by a coefficient vector: the value itself for a
# constant df, the sample average for a covariate-dependent df. One element
# for the Student-t copula, one per Student-t pair of a PCC.
.tc_mean_df <- function(par, object) {
  idx <- object$index
  V <- .tc_V(object)
  if (is.null(V)) return(2.05 + exp(par[idx$df]))
  dfl <- if (is.null(idx$dfl)) as.list(idx$df) else idx$dfl
  unlist(lapply(seq_along(dfl), function(k) {
    ii <- dfl[[k]]
    if (!length(ii)) return(NULL)
    if (.tc_df_const(V[[k]])) 2.05 + exp(par[[ii]]) else
      mean(2.05 + exp(drop(V[[k]] %*% par[ii])))
  }))
}

#' Dependence and scale parameters on their natural scale
#'
#' Delta-method standard errors and Wald intervals for the marginal
#' scale/shape parameters, the copula correlations (`rho12`, `rho13`,
#' `rho23`), the partial correlation `rho23|1` and, for the Student-t
#' copula, the degrees of freedom. Intervals for correlations are computed
#' on the Fisher-z scale and intervals for positive parameters on the log
#' scale, then back-transformed.
#'
#' @param object A `tricopula` fit.
#' @param level Confidence level.
#' @param type Covariance type, `"hessian"` or `"sandwich"`; see
#'   [tricopula-methods] for the default.
#' @return A data frame with estimates, standard errors and interval limits.
#' @export
tc_dependence <- function(object, level = 0.95, type = NULL) {
  .tc_check(object)
  type <- .tc_vcov_type(object, type)
  est <- .tc_natural(object$coefficients, object)
  se <- lo <- hi <- rep(NA_real_, length(est))
  if (!is.null(object$vcov)) {
    V <- stats::vcov(object, type = type)
    zq <- stats::qnorm(1 - (1 - level) / 2)
    J <- numDeriv::jacobian(.tc_natural, object$coefficients, object = object)
    na_par <- is.na(diag(V))
    V[is.na(V)] <- 0
    se <- sqrt(pmax(diag(J %*% V %*% t(J)), 0))
    is_rho <- grepl("^rho", names(est))
    is_df <- grepl("^df", names(est))
    if (any(na_par)) se[is_df] <- NA_real_
    # correlations: Fisher z; scale/shape: log; df: log(df - 2.05);
    # pair-copula parameters and Kendall's tau of a PCC: plain Wald
    lo[is_rho] <- tanh(atanh(est[is_rho]) - zq * se[is_rho] / (1 - est[is_rho]^2))
    hi[is_rho] <- tanh(atanh(est[is_rho]) + zq * se[is_rho] / (1 - est[is_rho]^2))
    plain <- grepl("^(par|tau)", names(est))
    lo[plain] <- est[plain] - zq * se[plain]
    hi[plain] <- est[plain] + zq * se[plain]
    pos <- !is_rho & !is_df & !plain
    lo[pos] <- est[pos] * exp(-zq * se[pos] / est[pos])
    hi[pos] <- est[pos] * exp(zq * se[pos] / est[pos])
    if (any(is_df)) {
      d <- est[is_df] - 2.05
      lo[is_df] <- 2.05 + d * exp(-zq * se[is_df] / d)
      hi[is_df] <- 2.05 + d * exp(zq * se[is_df] / d)
    }
  }
  data.frame(estimate = est, std.error = se, lower = lo, upper = hi,
             row.names = names(est))
}

#' @rdname tricopula-methods
#' @export
confint.tricopula <- function(object, parm, level = 0.95,
                              type = NULL, ...) {
  type <- .tc_vcov_type(object, type)
  b <- object$coefficients
  se <- sqrt(diag(stats::vcov(object, type = type)))
  zq <- stats::qnorm(1 - (1 - level) / 2)
  ci <- cbind(b - zq * se, b + zq * se)
  colnames(ci) <- paste0(format(100 * c((1 - level) / 2, 1 - (1 - level) / 2)), " %")
  if (!missing(parm)) ci <- ci[parm, , drop = FALSE]
  ci
}

#' @rdname tricopula-methods
#' @export
summary.tricopula <- function(object, type = NULL, ...) {
  type <- .tc_vcov_type(object, type)
  b <- object$coefficients
  se <- if (!is.null(object$vcov)) {
    sqrt(pmax(diag(stats::vcov(object, type = type)), 0))
  } else rep(NA_real_, length(b))
  z <- b / se
  tab <- cbind(Estimate = b, `Std. Error` = se, `z value` = z,
               `Pr(>|z|)` = 2 * stats::pnorm(abs(z), lower.tail = FALSE))
  structure(list(call = object$call, copula = object$copula,
    margins = object$margins, n = object$n, logLik = object$logLik,
    convergence = object$convergence, coefficients = tab,
    natural = tc_dependence(object, type = type), vcov_type = type,
    weights_type = object$weights_type,
    truncation = object$truncation, truncation_type = object$truncation_type,
    censored = object$censored,
    smooths = .tc_smooth_table(object), edf = if (!is.null(object$edf)) sum(object$edf),
    sp_converged = object$sp_converged,
    correlation = object$correlation, df = object$df,
    AIC = stats::AIC(object), BIC = stats::BIC(object),
    max_abs_gradient = max(abs(object$gradient))),
    class = "summary.tricopula")
}

#' @rdname tricopula-methods
#' @export
print.tricopula <- function(x, ...) {
  cat("TriCopula trivariate continuous regression\n")
  cat("Copula:", x$copula, "\n")
  cat("Margins:", paste(x$margins, collapse = ", "), "\n")
  cat("Observations:", x$n, "\n")
  .tc_print_ct(x)
  cat("Log-likelihood:", format(x$logLik), "\n")
  cat("Optimizer convergence code:", x$convergence, "\n")
  if (x$copula == "pcc") {
    pl <- .tc_pcc_labels(x$root)$short
    cat("Pair copulas (", paste0("c", pl, collapse = ", "), "): ",
        paste(x$families, collapse = ", "), " \n", sep = "")
    cat("Average Kendall's tau (", paste(pl, collapse = ", "), "): ",
        paste(format(colMeans(x$ktau), digits = 3), collapse = ", "), " \n", sep = "")
    if (isFALSE(x$simplified)) {
      cat("Non-simplified vine: c", pl[3], " depends on u", .tc_pcc_order(x$root)[1],
          " through a penalized spline (", x$nonsimplified$ncoef, " coefficients).\n", sep = "")
    }
    return(invisible(x))
  }
  if (x$copula == "student") {
    cat(if (isFALSE(x$df_constant)) "Average copula degrees of freedom:" else
      "Copula degrees of freedom:", x$df, "\n")
  }
  cat("Estimated copula correlation matrix:\n")
  print(x$correlation)
  invisible(x)
}

#' @rdname tricopula-methods
#' @export
print.summary.tricopula <- function(x, ...) {
  cat("TriCopula model summary\n")
  cat("Copula:", x$copula, " | Margins:", paste(x$margins, collapse = ", "),
      " | n:", x$n, "\n")
  cat("logLik:", format(x$logLik), " | AIC:", format(x$AIC),
      " | BIC:", format(x$BIC), "\n")
  cat("Convergence code:", x$convergence, " | max |gradient|:",
      format(x$max_abs_gradient, digits = 3), " | SE type:", x$vcov_type, "\n")
  if (!is.null(x$weights_type)) {
    cat("Weights:", x$weights_type,
        if (x$weights_type == "sampling") "(sandwich standard errors by default)" else
          "(rows count as replicates; model-based standard errors by default)", "\n")
  }
  .tc_print_ct(x)
  cat("\n")
  cat("Coefficients (unconstrained scale):\n")
  stats::printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE)
  cat("\nScale, shape and dependence parameters (natural scale, delta method):\n")
  print(round(x$natural, 4))
  if (!is.null(x$smooths)) {
    cat("\nPenalized smooth terms (total effective degrees of freedom ",
        format(x$edf, digits = 4), "; standard errors above are Bayesian;\n",
        "the tests are approximate):\n", sep = "")
    sm <- x$smooths
    num <- vapply(sm, is.numeric, logical(1))
    sm[num] <- lapply(sm[num], signif, 4)
    print(sm, row.names = FALSE)
    if (isFALSE(x$sp_converged)) cat("The smoothing-parameter iteration did not converge.\n")
  }
  invisible(x)
}
