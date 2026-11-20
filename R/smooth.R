# Penalized smooth terms (0.12.0, generalized in 0.12.1) --------------------
#
# A formula of any predictor may contain mgcv-style smooths: s(), te(), ti(),
# with or without a `by` variable. The basis and the penalties of each smooth
# come from mgcv::smoothCon() with the identifiability constraint absorbed.
# The model is fitted by maximizing
#     l(theta) - 1/2 sum_j lambda_j theta' S_j theta,
# where a smooth may have several penalties S_j on the same coefficients
# (tensor products, adaptive or shrinkage bases). The smoothing parameters
# are chosen by the generalized Fellner-Schall update of Wood and Fasiolo
# (2017, Biometrics 73, 1071-1081):
#     lambda_j <- lambda_j * (tr(S_lambda^- S_j) - tr(V_p S_j)) /
#                 (theta' S_j theta),   V_p = (H + S_lambda)^-1,
# with H the observed information of the unpenalized log-likelihood and
# S_lambda^- the pseudo-inverse of the total penalty. For a smooth with one
# penalty tr(S_lambda^- S_j) = rank(S_j) / lambda_j.
# Formulas without smooth terms never reach this file.

.tc_smooth_funs <- c("s", "te", "ti", "t2")

.tc_has_smooth <- function(formula) {
  any(.tc_smooth_funs %in% all.names(formula))
}

# Split a formula into its parametric part and its smooths and build the
# smooth bases on `data`. Returns list(pf, smooths); each smooth is the mgcv
# smooth object (needed by PredictMat) with X (columns named after the
# smooth's label) and tc_S, the list of its penalty matrices.
# `basis`: smooth objects of an earlier fit of the same formula; their
# knots and constraints are then reused and only X is evaluated on `data`.
.tc_smooth_setup <- function(formula, data, basis = NULL) {
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("Package 'mgcv' is required for formulas with smooth terms.", call. = FALSE)
  }
  if ("t2" %in% all.names(formula)) {
    stop("t2() smooths are not supported; use te() or ti().", call. = FALSE)
  }
  one_sided <- length(formula) == 2L
  f <- if (one_sided) {
    stats::as.formula(paste(".tc_dummy ~", paste(deparse(formula[[2]]), collapse = " ")),
                      env = environment(formula))
  } else formula
  gp <- mgcv::interpret.gam(f)
  pf <- gp$pf
  if (one_sided) pf <- pf[-2]
  environment(pf) <- environment(formula)
  if (!is.null(basis)) {
    smooths <- lapply(basis, function(sm) {
      sm$X <- mgcv::PredictMat(sm, data)
      colnames(sm$X) <- sm$tc_names
      sm
    })
    return(list(pf = pf, smooths = smooths))
  }
  smooths <- list()
  for (spec in gp$smooth.spec) {
    need <- c(spec$term, if (spec$by != "NA") spec$by)
    if (!all(need %in% names(data))) {
      stop("Variables of ", spec$label, " must be columns of the data.", call. = FALSE)
    }
    # a factor `by` variable gives one smooth per level
    sl <- mgcv::smoothCon(spec, data = data, knots = NULL, absorb.cons = TRUE,
                          scale.penalty = TRUE)
    for (sm in sl) {
      sm$tc_names <- paste0(sm$label, ".", seq_len(ncol(sm$X)))
      colnames(sm$X) <- sm$tc_names
      sm$tc_S <- lapply(sm$S, function(S) (S + t(S)) / 2)
      smooths[[length(smooths) + 1L]] <- sm
    }
  }
  list(pf = pf, smooths = smooths)
}

# Smooth columns of a design matrix for new data.
.tc_smooth_predict <- function(smooths, newdata) {
  do.call(cbind, lapply(smooths, function(sm) {
    P <- mgcv::PredictMat(sm, newdata)
    colnames(P) <- sm$tc_names
    P
  }))
}

# Append the smooth columns to a parametric design matrix; records in each
# smooth the positions of its columns inside the predictor (tc_cols).
.tc_smooth_bind <- function(X, smooths) {
  pos <- ncol(X)
  for (i in seq_along(smooths)) {
    k <- ncol(smooths[[i]]$X)
    smooths[[i]]$tc_cols <- pos + seq_len(k)
    pos <- pos + k
    X <- cbind(X, smooths[[i]]$X)
    smooths[[i]]$X <- NULL                 # kept in the design matrix only
  }
  list(X = X, smooths = smooths)
}

# Collect the smooths of all predictors. `blocks` is a list of
# list(prefix, positions (of the predictor's coefficients in the parameter
# vector), smooths). Returns
#   sblocks: one entry per smooth: label, ind (its coefficients), smooth
#   pens:    one entry per penalty: ind, S, block (index into sblocks), label
.tc_collect_penalties <- function(blocks) {
  sblocks <- list(); pens <- list()
  for (b in blocks) {
    for (sm in b$smooths) {
      lab <- paste0(b$prefix, ":", sm$label)
      ind <- b$positions[sm$tc_cols]
      sblocks[[length(sblocks) + 1L]] <- list(label = lab, ind = ind, smooth = sm)
      m <- length(sm$tc_S)
      for (j in seq_len(m)) {
        pens[[length(pens) + 1L]] <- list(
          ind = ind, S = sm$tc_S[[j]], block = length(sblocks),
          label = if (m == 1L) lab else paste0(lab, ".", j))
      }
    }
  }
  list(sblocks = sblocks, pens = pens)
}

.tc_pen_value <- function(par, pens, lam) {
  v <- 0
  for (k in seq_along(pens)) {
    b <- par[pens[[k]]$ind]
    v <- v + 0.5 * lam[k] * sum(b * (pens[[k]]$S %*% b))
  }
  v
}
.tc_pen_grad <- function(par, pens, lam) {
  g <- numeric(length(par))
  for (k in seq_along(pens)) {
    ii <- pens[[k]]$ind
    g[ii] <- g[ii] + lam[k] * drop(pens[[k]]$S %*% par[ii])
  }
  g
}
.tc_pen_matrix <- function(npar, pens, lam) {
  S <- matrix(0, npar, npar)
  for (k in seq_along(pens)) {
    ii <- pens[[k]]$ind
    S[ii, ii] <- S[ii, ii] + lam[k] * pens[[k]]$S
  }
  S
}

# tr(S_lambda^- S_j) for every penalty: within each smooth the pseudo-inverse
# of its total penalty sum_j lambda_j S_j.
.tc_pen_trace <- function(pens, lam) {
  tr <- numeric(length(pens))
  blk <- vapply(pens, `[[`, integer(1), "block")
  for (b in unique(blk)) {
    kk <- which(blk == b)
    if (length(kk) == 1L) {
      ev <- eigen(pens[[kk]]$S, symmetric = TRUE, only.values = TRUE)$values
      tr[kk] <- sum(ev > max(ev) * 1e-8) / lam[kk]
      next
    }
    St <- Reduce(`+`, lapply(kk, function(k) lam[k] * pens[[k]]$S))
    e <- eigen(St, symmetric = TRUE)
    keep <- e$values > max(e$values) * 1e-8
    Sinv <- e$vectors[, keep, drop = FALSE] %*%
      (t(e$vectors[, keep, drop = FALSE]) / e$values[keep])
    for (k in kk) tr[k] <- sum(Sinv * pens[[k]]$S)
  }
  tr
}

# One generalized Fellner-Schall update of the smoothing parameters. H is
# the unpenalized observed information at the current penalized estimate.
# The step is limited to a factor of 20 either way and lambda is kept in
# [1e-7, 1e7]; NULL is returned if H + S_lambda cannot be inverted. Also
# returns the effective degrees of freedom of each smooth at the current
# lambda, which the caller uses as its stopping rule: the smoothing
# parameter of an effect that is really linear runs to infinity slowly,
# while its effective degrees of freedom have long settled.
.tc_fs_update <- function(par, H, pens, lam, nblocks) {
  Hp <- H + .tc_pen_matrix(length(par), pens, lam)
  Vp <- try(solve(Hp), silent = TRUE)
  if (inherits(Vp, "try-error") || any(!is.finite(Vp))) return(NULL)
  trS <- .tc_pen_trace(pens, lam)
  new <- lam
  edf <- numeric(nblocks)
  for (k in seq_along(pens)) {
    ii <- pens[[k]]$ind
    b <- par[ii]
    den <- sum(b * (pens[[k]]$S %*% b))
    num <- trS[k] - sum(Vp[ii, ii] * pens[[k]]$S)
    fac <- if (!is.finite(num) || num <= 0) 1 / 20 else if (den <= 0) 20 else num / den
    new[k] <- lam[k] * min(max(fac, 1 / 20), 20)
    # effective degrees of freedom of the smooth at the current lambda
    edf[pens[[k]]$block] <- sum(Vp[ii, , drop = FALSE] * H[ii, , drop = FALSE])
  }
  list(lam = pmin(pmax(new, 1e-7), 1e7), edf = edf)
}

# Newton iteration with a fixed matrix M in place of the Hessian (the
# penalized Hessian built from the information of the previous smoothing
# step): each step costs one gradient, no numerical Hessian. Steps are
# accepted only if they do not increase the objective. ok = FALSE if the
# gradient is not below `tol` after `maxit` steps.
.tc_newton_fixed <- function(par, nll, ngr, M, maxit = 12L, tol = 1e-5) {
  R <- try(chol(M), silent = TRUE)
  if (inherits(R, "try-error")) return(list(par = par, ok = FALSE))
  val <- nll(par)
  for (k in seq_len(maxit)) {
    g <- ngr(par)
    if (max(abs(g)) < tol) return(list(par = par, ok = TRUE))
    step <- backsolve(R, backsolve(R, g, transpose = TRUE))
    lamb <- 1; improved <- FALSE
    for (h in 1:8) {
      cand <- par - lamb * step
      cval <- nll(cand)
      if (is.finite(cval) && cval <= val) {
        par <- cand; val <- cval; improved <- TRUE
        break
      }
      lamb <- lamb / 2
    }
    if (!improved) break
  }
  list(par = par, ok = max(abs(ngr(par))) < tol)
}

# Smooths of a fitted model: one row per smooth (used by summary() and
# tc_smooth()). `lambda` is the smoothing parameter of the first penalty of
# the smooth; all of them are in object$sp.
.tc_smooth_table <- function(object, V = object$vcov) {
  sb <- object$smooth_blocks
  if (is.null(sb)) return(NULL)
  b <- object$coefficients
  blk <- vapply(object$penalties, `[[`, integer(1), "block")
  out <- data.frame(
    term = vapply(sb, `[[`, "", "label"),
    coefficients = vapply(sb, function(p) length(p$ind), integer(1)),
    penalties = vapply(seq_along(sb), function(i) sum(blk == i), integer(1)),
    lambda = vapply(seq_along(sb), function(i) {
      k <- which(blk == i)
      if (length(k)) unname(object$sp[k[1]]) else NA_real_
    }, numeric(1)),
    edf = vapply(sb, function(p) {
      if (is.null(object$edf)) NA_real_ else sum(object$edf[p$ind])
    }, numeric(1)),
    stringsAsFactors = FALSE)
  # approximate Wald test of the whole smooth: edf (not less than 1) as
  # degrees of freedom
  if (!is.null(V)) {
    st <- vapply(sb, function(p) {
      r <- try(sum(b[p$ind] * solve(V[p$ind, p$ind, drop = FALSE], b[p$ind])),
               silent = TRUE)
      if (inherits(r, "try-error")) NA_real_ else r
    }, numeric(1))
    out$chi.sq <- st
    out$p.value <- stats::pchisq(st, pmax(out$edf, 1), lower.tail = FALSE)
  }
  rownames(out) <- NULL
  out
}

# Covariance matrix corrected for the uncertainty of the smoothing
# parameters (Wood, Pya and Saefken, 2016, JASA 111, 1548-1563, first-order
# correction): V_c = V_p + J V_rho J', with J = d theta-hat / d rho,
# rho = log lambda, column j equal to -lambda_j V_p S_j theta-hat, and V_rho
# the inverse of the negative Hessian of the Laplace-approximate marginal
# likelihood in rho. That Hessian is obtained by forward differences of the
# gradient
#   g_j = lambda_j / 2 * (tr(S_lambda^- S_j) - tr(V_p S_j) - theta' S_j theta)
# (the same approximation as the Fellner-Schall update: the dependence of H
# on rho is neglected), which takes one refit per smoothing parameter.
# Directions in which the criterion is flat (smoothing parameters at their
# bounds, effects shrunk to the null space) are left out.
.tc_sp_correction <- function(object, h = 0.2) {
  pens <- object$penalties
  if (is.null(pens)) stop("The fit has no penalized smooth terms.")
  if (is.null(object$vcov) || anyNA(object$vcov)) {
    stop("The fit has no complete covariance matrix.")
  }
  npar <- length(object$coefficients)
  grad_rho <- function(par, lam, Hu) {
    Vp <- solve((Hu + t(Hu)) / 2 + .tc_pen_matrix(npar, pens, lam))
    trS <- .tc_pen_trace(pens, lam)
    vapply(seq_along(pens), function(k) {
      ii <- pens[[k]]$ind
      b <- par[ii]
      0.5 * lam[k] * (trS[k] - sum(Vp[ii, ii] * pens[[k]]$S) - sum(b * (pens[[k]]$S %*% b)))
    }, numeric(1))
  }
  hess_u <- function(fit) {
    -numDeriv::jacobian(function(p) tc_gradient(fit, p), fit$coefficients,
                        method = "simple")
  }
  lam <- unname(object$sp)
  par <- unname(object$coefficients)
  Vp <- object$vcov
  g0 <- grad_rho(par, lam, hess_u(object))
  m <- length(pens)
  A <- matrix(0, m, m)
  J <- matrix(0, npar, m)
  for (k in seq_len(m)) {
    ii <- pens[[k]]$ind
    J[, k] <- -lam[k] * drop(Vp[, ii, drop = FALSE] %*% (pens[[k]]$S %*% par[ii]))
    l1 <- lam
    l1[k] <- lam[k] * exp(h)
    fk <- .tc_refit(object, object$data, sp = l1, hessian = FALSE)
    if (inherits(fk, "try-error")) stop("A refit for the smoothing-parameter correction failed.")
    A[, k] <- (grad_rho(unname(fk$coefficients), l1, hess_u(fk)) - g0) / h
  }
  M <- -(A + t(A)) / 2
  e <- eigen(M, symmetric = TRUE)
  keep <- e$values > max(max(e$values) * 1e-4, 1e-8)
  Vrho <- if (any(keep)) {
    e$vectors[, keep, drop = FALSE] %*% (t(e$vectors[, keep, drop = FALSE]) / e$values[keep])
  } else matrix(0, m, m)
  Vc <- Vp + J %*% Vrho %*% t(J)
  dimnames(Vc) <- dimnames(Vp)
  attr(Vc, "V_rho") <- Vrho
  Vc
}

#' Fitted smooth terms with pointwise intervals
#'
#' Evaluates one penalized smooth of a fitted model on a grid, on the scale
#' of its linear predictor, with pointwise intervals. The smooth is centred
#' (it sums to zero over the estimation sample); the intercept of the
#' predictor is not included.
#'
#' @param object A `tricopula` fit with at least one smooth term.
#' @param which Number or label of the smooth, as listed in
#'   `summary(object)$smooths`.
#' @param n Number of grid points over the observed range of the covariate.
#' @param newdata Optional data frame with the covariates of the smooth; if
#'   given it replaces the grid. Required for smooths of several covariates
#'   and for smooths with a `by` variable.
#' @param level Confidence level of the pointwise intervals.
#' @param type `"bayesian"` (the default) uses the Bayesian covariance
#'   matrix, which is conditional on the selected smoothing parameters;
#'   `"corrected"` adds the first-order correction for their uncertainty
#'   (see [tricopula-methods]; it refits the model once per smoothing
#'   parameter).
#' @param vcov Optional covariance matrix to use instead, e.g. the result of
#'   `vcov(object, type = "corrected")` computed once.
#' @return A data frame with the covariate(s), `fit`, `se`, `lower`,
#'   `upper`; attribute `"term"` holds the label of the smooth.
#' @export
tc_smooth <- function(object, which = 1L, n = 100L, newdata = NULL, level = 0.95,
                      type = c("bayesian", "corrected"), vcov = NULL) {
  .tc_check(object)
  type <- match.arg(type)
  sb <- object$smooth_blocks
  if (is.null(sb)) stop("The fit has no s() terms.")
  labs <- vapply(sb, `[[`, "", "label")
  k <- if (is.character(which)) match(which, labs) else as.integer(which)
  if (length(k) != 1L || is.na(k) || k < 1L || k > length(sb)) {
    stop("which must be one of: ", paste(labs, collapse = ", "), " (or its number).")
  }
  sm <- sb[[k]]$smooth
  if (isTRUE(sm$tc_ns)) {
    # spline in the conditioning value of a non-simplified vine
    if (is.null(newdata)) newdata <- data.frame(u = seq(0.01, 0.99, length.out = n))
    if (!is.numeric(newdata$u) || any(newdata$u < 0 | newdata$u > 1)) {
      stop("Supply newdata with a column u of values in [0, 1] for this smooth.")
    }
    P <- object$nonsimplified$basis(newdata$u, FALSE)$B
    vars <- "u"
  } else {
  vars <- c(sm$term, if (sm$by != "NA") sm$by)
  if (is.null(newdata)) {
    if (length(vars) != 1L) {
      stop("Supply newdata with the columns ", paste(vars, collapse = ", "),
           " for this smooth.")
    }
    x <- object$data[[sm$term]]
    newdata <- stats::setNames(data.frame(seq(min(x), max(x), length.out = n)), sm$term)
  }
  P <- mgcv::PredictMat(sm, newdata)
  }
  ii <- sb[[k]]$ind
  fit <- drop(P %*% object$coefficients[ii])
  V <- if (!is.null(vcov)) vcov else if (type == "corrected") {
    stats::vcov(object, type = "corrected")
  } else object$vcov
  se <- if (is.null(V)) rep(NA_real_, length(fit)) else {
    sqrt(pmax(rowSums((P %*% V[ii, ii, drop = FALSE]) * P), 0))
  }
  zq <- stats::qnorm(1 - (1 - level) / 2)
  out <- cbind(newdata[vars], fit = fit, se = se, lower = fit - zq * se,
               upper = fit + zq * se)
  attr(out, "term") <- labs[k]
  out
}
