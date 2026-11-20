# Post-estimation ---------------------------------------------------------

# Design matrix of a predictor for new data: an intercept column if the
# parametric part has no covariates, plus the columns of any s() terms.
.tc_add_smooth <- function(Xn, smooths, newdata) {
  if (nrow(Xn) != nrow(newdata)) Xn <- matrix(1, nrow(newdata), 1)
  if (is.null(smooths)) Xn else cbind(Xn, .tc_smooth_predict(smooths, newdata))
}

.tc_new_eta <- function(object, newdata) {
  idx <- object$index
  lapply(1:3, function(j) {
    Xn <- if (is.null(newdata)) object$design[[j]] else {
      mf <- stats::model.frame(object$terms[[j]], newdata,
                               na.action = stats::na.fail,
                               xlev = object$xlevels[[j]])
      .tc_add_smooth(stats::model.matrix(object$terms[[j]], mf,
                                         contrasts.arg = object$contrasts[[j]]),
                     object$smooths$mean[[j]], newdata)
    }
    drop(Xn %*% object$coefficients[idx$beta[[j]]])
  })
}

# Log scale / log shape of each margin for each covariate row.
.tc_new_extra <- function(object, newdata) {
  idx <- object$index
  lapply(1:3, function(j) {
    Zn <- if (is.null(newdata)) object$scale_design[[j]] else {
      mf <- stats::model.frame(object$scale_terms[[j]], newdata,
                               na.action = stats::na.fail,
                               xlev = object$scale_xlevels[[j]])
      Zm <- stats::model.matrix(object$scale_terms[[j]], mf,
                                contrasts.arg = object$scale_contrasts[[j]])
      .tc_add_smooth(Zm, object$smooths$scale[[j]], newdata)
    }
    drop(Zn %*% object$coefficients[idx$extral[[j]]])
  })
}

# Third (shape) predictor of each three-parameter margin for each covariate
# row; NULL for two-parameter margins.
.tc_new_shape <- function(object, newdata) {
  idx <- object$index
  lapply(1:3, function(j) {
    ii <- idx$shapel[[j]]
    if (!length(ii)) return(NULL)
    Qn <- if (is.null(newdata)) object$shape_design[[j]] else {
      mf <- stats::model.frame(object$shape_terms[[j]], newdata,
                               na.action = stats::na.fail,
                               xlev = object$shape_xlevels[[j]])
      Qm <- stats::model.matrix(object$shape_terms[[j]], mf,
                                contrasts.arg = object$shape_contrasts[[j]])
      .tc_add_smooth(Qm, object$smooths$shape[[j]], newdata)
    }
    drop(Qn %*% object$coefficients[ii])
  })
}

# Dependence predictors (list of three vectors) for each covariate row.
.tc_new_dep_eta <- function(object, newdata) {
  idx <- object$index
  lapply(1:3, function(k) {
    Wn <- if (is.null(newdata)) object$dep_design[[k]] else {
      mf <- stats::model.frame(object$dep_terms[[k]], newdata,
                               na.action = stats::na.fail,
                               xlev = object$dep_xlevels[[k]])
      Wm <- stats::model.matrix(object$dep_terms[[k]], mf,
                                contrasts.arg = object$dep_contrasts[[k]])
      .tc_add_smooth(Wm, object$smooths$dep[[k]], newdata)
    }
    drop(Wn %*% object$coefficients[idx$depl[[k]]])
  })
}

# log(df - 2.05) of each pair of a PCC (NULL for non-Student pairs).
.tc_pcc_tau <- function(object, newdata = NULL) {
  if (!isFALSE(object$df_constant)) {
    it <- .tc_pair_base(object$families) == "student"
    pos <- cumsum(it)
    return(lapply(1:3, function(k) {
      if (it[k]) object$coefficients[[object$index$df[pos[k]]]] else NULL
    }))
  }
  .tc_new_tau(object, newdata)
}

# log(df - 2.05) predictors of a fit with covariate-dependent df: a list
# with one element for the Student-t copula and three for a PCC (NULL for
# pairs without df; a scalar for a pair whose df is constant).
.tc_new_tau <- function(object, newdata = NULL) {
  idx <- object$index
  lapply(seq_along(idx$dfl), function(k) {
    ii <- idx$dfl[[k]]
    if (!length(ii)) return(NULL)
    Vk <- object$df_design[[k]]
    if (.tc_df_const(Vk)) return(object$coefficients[[ii]])
    if (!is.null(newdata)) {
      mf <- stats::model.frame(object$df_terms[[k]], newdata,
                               na.action = stats::na.fail,
                               xlev = object$df_xlevels[[k]])
      Vk <- .tc_add_smooth(stats::model.matrix(object$df_terms[[k]], mf,
                                               contrasts.arg = object$df_contrasts[[k]]),
                           object$smooths$df[[k]], newdata)
    }
    drop(Vk %*% object$coefficients[ii])
  })
}

# Fitted correlations (rho12, rho13, rho23, rho23|1) for each covariate row.
.tc_new_rho <- function(object, newdata) {
  if (is.null(newdata)) return(object$rho)
  eta <- .tc_new_dep_eta(object, newdata)
  dep <- .tc_dep_vec(eta[[1]], eta[[2]], eta[[3]])
  cbind(rho12 = dep$r12, rho13 = dep$r13, rho23 = dep$r23, `rho23|1` = dep$p23)
}

#' Predict marginal distribution parameters, means or dependence
#'
#' @param object A `tricopula` fit.
#' @param newdata Optional data frame; defaults to the estimation sample.
#' @param type `"parameters"` returns the marginal distribution parameters;
#'   `"mean"` returns the conditional means `E(Y_j | x)`; `"dependence"`
#'   returns the fitted copula correlations, the partial correlation
#'   `rho23|1` and the implied Kendall's tau for each row.
#' @param ... Unused.
#' @return For `type = "parameters"`, a list of three data frames; for
#'   `type = "mean"`, a data frame with one column per response; for
#'   `type = "dependence"`, a data frame with columns `rho12`, `rho13`,
#'   `rho23`, `rho23|1`, `tau12`, `tau13`, `tau23`. For a PCC the columns
#'   are the pair parameters and Kendall's tau of the three pairs, named
#'   after the pairs of the vine (e.g. `par13|2` with `root = 2`). If the
#'   degrees of freedom depend on covariates (`df` formula), columns with
#'   the fitted degrees of freedom are appended.
#' @export
predict.tricopula <- function(object, newdata = NULL,
                              type = c("parameters", "mean", "dependence"), ...) {
  type <- match.arg(type)
  if (type == "dependence" && object$copula == "pcc") {
    eta <- .tc_new_dep_eta(object, newdata)
    fam <- object$families
    pl <- .tc_pcc_labels(object$root)$short
    out <- cbind(sapply(1:3, function(k) .tc_pair_natural(eta[[k]], fam[k])),
                 sapply(1:3, function(k) .tc_pair_ktau(eta[[k]], fam[k])))
    out <- matrix(out, ncol = 6, dimnames = list(NULL, c(
      paste0("par", pl), paste0("tau", pl))))
    # non-simplified vine: fitted values of the conditional pair for the
    # estimation sample; for new data its value averaged over the
    # conditioning value on the predictor scale
    if (!is.null(object$nonsimplified) && is.null(newdata)) {
      out[, 3] <- object$pair[, 3]; out[, 6] <- object$ktau[, 3]
    }
    if (isFALSE(object$df_constant)) {
      tl <- .tc_new_tau(object, newdata)
      for (k in which(lengths(tl) > 0L)) {
        out <- cbind(out, rep_len(2.05 + exp(tl[[k]]), nrow(out)))
        colnames(out)[ncol(out)] <- paste0("df", pl[k])
      }
    }
    return(as.data.frame(out, optional = TRUE))
  }
  if (type == "dependence") {
    rho <- .tc_new_rho(object, newdata)
    tau <- 2 / pi * asin(rho[, 1:3, drop = FALSE])
    colnames(tau) <- c("tau12", "tau13", "tau23")
    out <- cbind(rho, tau)
    if (isFALSE(object$df_constant)) {
      out <- cbind(out, df = rep_len(2.05 + exp(.tc_new_tau(object, newdata)[[1]]),
                                     nrow(out)))
    }
    return(as.data.frame(out, optional = TRUE))
  }
  eta <- .tc_new_eta(object, newdata)
  extra <- .tc_new_extra(object, newdata)
  shp <- .tc_new_shape(object, newdata)
  if (type == "mean") {
    tb <- if (identical(object$truncation_type, "margin")) {
      .tc_new_trunc(object, newdata, length(eta[[1]]))
    }
    out <- lapply(1:3, function(j) {
      if (is.null(tb[[j]])) {
        return(.tc_margin_mean(eta[[j]], extra[[j]], object$margins[j], shp[[j]]))
      }
      # mean of the truncated margin: integral of its quantile function
      Fab <- .tc_trunc_F(tb[[j]], eta[[j]], extra[[j]], object$margins[j], shp[[j]])
      gl <- .tc_gl(64L)
      G <- gl$x^3 * (10 - 15 * gl$x + 6 * gl$x^2)
      gw <- gl$w * 30 * gl$x^2 * (1 - gl$x)^2
      m <- numeric(length(eta[[j]]))
      for (q in seq_along(G)) {
        m <- m + gw[q] * .tc_margin_quantile(.tc_clip(Fab$Fa + (Fab$Fb - Fab$Fa) * G[q]),
                                             eta[[j]], extra[[j]], object$margins[j], shp[[j]])
      }
      m
    })
    names(out) <- object$response_names
    return(as.data.frame(out, optional = TRUE))
  }
  out <- lapply(1:3, function(j) {
    .tc_margin_pars(eta[[j]], extra[[j]], object$margins[j], shp[[j]])
  })
  names(out) <- paste0("margin", 1:3)
  out
}

#' Simulate responses from a fitted model
#'
#' @param object A `tricopula` fit.
#' @param nsim Number of simulated data sets.
#' @param seed Optional random seed.
#' @param newdata Optional data frame of covariates; defaults to the
#'   estimation sample.
#' @param ... Unused.
#' @return A data frame of simulated responses (a list of data frames if
#'   `nsim > 1`).
#' @export
simulate.tricopula <- function(object, nsim = 1, seed = NULL, newdata = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)
  eta <- .tc_new_eta(object, newdata)
  extra <- .tc_new_extra(object, newdata)
  shp <- .tc_new_shape(object, newdata)
  n <- length(eta[[1]])
  pcc <- object$copula == "pcc"
  rho <- if (pcc) NULL else .tc_new_rho(object, newdata)
  dep_eta <- if (pcc) .tc_new_dep_eta(object, newdata)
  ptau <- if (pcc) .tc_pcc_tau(object, newdata)
  ord <- if (pcc) .tc_pcc_order(object$root)
  dfc <- !isFALSE(object$df_constant) || object$copula != "student"
  dfv <- if (dfc) object$df else 2.05 + exp(.tc_new_tau(object, newdata)[[1]])
  if (!pcc && !dfc && nrow(rho) != n) rho <- rho[rep_len(1L, n), , drop = FALSE]
  one <- function() {
    U <- if (pcc) {
      Uv <- .tc_rpcc(n, dep_eta, object$families, ptau, object$nonsimplified,
                     object$coefficients[object$index$ns])
      if (!identical(ord, 1:3)) Uv[, ord] <- Uv
      Uv
    } else if (object$dep_constant && dfc) {
      tc_rcopula(n, object$correlation, object$copula, object$df)
    } else {
      .tc_rcopula_vec(n, rho[, 1], rho[, 2], rho[, 4], object$copula, dfv)
    }
    out <- lapply(1:3, function(j) {
      p <- U[, j]
      if (!is.null(tb[[j]]) && !joint) {
        Fab <- .tc_trunc_F(tb[[j]], eta[[j]], extra[[j]], object$margins[j], shp[[j]])
        p <- .tc_clip(Fab$Fa + (Fab$Fb - Fab$Fa) * p)
      }
      .tc_margin_quantile(p, eta[[j]], extra[[j]], object$margins[j], shp[[j]])
    })
    names(out) <- object$response_names
    as.data.frame(out, optional = TRUE)
  }
  # truncation: type "margin" draws from the truncated margins; type "joint"
  # redraws the rows whose truncated response falls outside its bounds
  tb <- .tc_new_trunc(object, newdata, n)
  joint <- identical(object$truncation_type, "joint")
  draw <- one
  if (joint) {
    jt <- which(!vapply(tb, is.null, logical(1)))
    draw <- function() {
      out <- one()
      for (it in 1:500) {
        bad <- out[[jt]] < tb[[jt]]$lower | out[[jt]] > tb[[jt]]$upper
        if (!any(bad)) break
        out[bad, ] <- one()[bad, ]
      }
      if (any(bad)) warning("Some simulated rows are outside the truncation bounds.", call. = FALSE)
      out
    }
  }
  if (nsim == 1) draw() else replicate(nsim, draw(), simplify = FALSE)
}

#' Probability integral transforms and per-observation log scores
#'
#' @param object A `tricopula` fit.
#' @return `tc_pit()` returns a data frame of the three fitted PITs
#'   `F_j(y_j | x)`; `tc_joint_logscore()` the per-observation joint
#'   log-density.
#' @export
tc_pit <- function(object) {
  .tc_check(object)
  as.data.frame(object$uniforms)
}

#' @rdname tc_pit
#' @export
tc_joint_logscore <- function(object) {
  .tc_check(object)
  object$loglik_i
}

#' Kendall's tau and tail dependence implied by the fitted copula
#'
#' With covariate-dependent dependence both functions use the
#' sample-average correlation matrix (and the sample-average degrees of
#' freedom if these depend on covariates); use `predict(type = "dependence")`
#' for observation-specific values. For a PCC, `tc_taildep()` reports the two
#' unconditional pairs of the vine (the pairs that contain the root). The
#' Frank copula has no tail dependence; the Joe copula has upper-tail
#' dependence `2 - 2^(1 / theta)`.
#' For elliptical copulas `tau = (2 / pi) asin(rho)`. The Gaussian copula has
#' zero tail dependence; for the Student-t copula the (equal) upper and
#' lower coefficients are
#' `2 * pt(-sqrt((df + 1) * (1 - rho) / (1 + rho)), df + 1)`.
#'
#' @param object A `tricopula` fit.
#' @return A 3 by 3 matrix.
#' @export
tc_tau <- function(object) {
  .tc_check(object)
  if (object$copula == "pcc") return(colMeans(object$ktau))
  tau <- 2 / pi * asin(object$correlation)
  dimnames(tau) <- list(object$response_names, object$response_names)
  tau
}

#' @rdname tc_tau
#' @export
tc_taildep <- function(object) {
  .tc_check(object)
  if (object$copula == "pcc") {
    eta <- vapply(1:3, function(k) {
      mean(drop(object$dep_design[[k]] %*%
                  object$coefficients[object$index$depl[[k]]]))
    }, numeric(1))
    ptau <- lapply(.tc_pcc_tau(object), function(t) if (is.null(t)) NULL else mean(t))
    td <- t(sapply(1:2, function(k) {
      .tc_pair_taildep(eta[k], object$families[k], ptau[[k]])
    }))
    rn <- object$response_names[.tc_pcc_order(object$root)]
    return(data.frame(pair = c(paste(rn[1], rn[2], sep = "-"), paste(rn[1], rn[3], sep = "-")),
                      family = object$families[1:2], lower = td[, "lower"],
                      upper = td[, "upper"]))
  }
  R <- object$correlation
  lam <- if (object$copula == "gaussian") matrix(0, 3, 3) else
    2 * stats::pt(-sqrt((object$df + 1) * (1 - R) / (1 + R)), object$df + 1)
  diag(lam) <- 1
  dimnames(lam) <- list(object$response_names, object$response_names)
  lam
}

.tc_gauss_legendre <- function(m) {
  i <- seq_len(m - 1L)
  J <- matrix(0, m, m)
  b <- i / sqrt(4 * i^2 - 1)
  J[cbind(i, i + 1L)] <- b
  J[cbind(i + 1L, i)] <- b
  e <- eigen(J, symmetric = TRUE)
  list(x = (e$values + 1) / 2, w = e$vectors[1, ]^2)
}

# Trivariate copula CDF at one point u. The Student-t case integrates the
# Gaussian orthant probability over the chi-square mixing variable, which
# allows non-integer degrees of freedom.
.tc_copula_cdf <- function(u, R, df, copula, gl = .tc_gauss_legendre(48L)) {
  pn <- function(upper) {
    as.numeric(mvtnorm::pmvnorm(lower = rep(-Inf, 3), upper = upper, corr = R,
                                algorithm = mvtnorm::TVPACK(abseps = 1e-8)))
  }
  if (copula == "gaussian") return(pn(stats::qnorm(u)))
  z <- stats::qt(u, df = df)
  sc <- sqrt(stats::qchisq(gl$x, df = df) / df)
  sum(gl$w * vapply(sc, function(s) pn(z * s), numeric(1)))
}

#' Joint probabilities from the fitted trivariate distribution
#'
#' Computes `P(Y1 <= t1, Y2 <= t2, Y3 <= t3 | x)` for each covariate row.
#' Requires the \pkg{mvtnorm} package.
#'
#' @param object A `tricopula` fit.
#' @param thresholds Numeric vector of three thresholds `(t1, t2, t3)`.
#' @param newdata Optional data frame; defaults to the estimation sample.
#' @return A numeric vector of joint probabilities.
#' @export
tc_jointprob <- function(object, thresholds, newdata = NULL) {
  .tc_check(object)
  if (!requireNamespace("mvtnorm", quietly = TRUE)) {
    stop("Package 'mvtnorm' is required for tc_jointprob().")
  }
  if (!is.numeric(thresholds) || length(thresholds) != 3L) {
    stop("thresholds must be a numeric vector of length 3.")
  }
  eta <- .tc_new_eta(object, newdata)
  extra <- .tc_new_extra(object, newdata)
  shp <- .tc_new_shape(object, newdata)
  U <- sapply(1:3, function(j) {
    .tc_clip(.tc_margin_cdf(thresholds[j], eta[[j]], extra[[j]], object$margins[j],
                            shp[[j]]))
  })
  U <- matrix(U, ncol = 3)
  # truncation (0.16.0). Type "margin": the PIT of the truncated margin.
  # Type "joint": P(Y <= t | a < Y_j < b) = (C(.., min(F_j(t_j), F_j(b)), ..) -
  # C(.., F_j(a), ..)) / (F_j(b) - F_j(a)), computed by two calls.
  tb <- .tc_new_trunc(object, newdata, nrow(U))
  if (!is.null(tb)) {
    jt <- which(!vapply(tb, is.null, logical(1)))
    Fab <- lapply(jt, function(j) {
      .tc_trunc_F(tb[[j]], eta[[j]], extra[[j]], object$margins[j], shp[[j]])
    })
    if (identical(object$truncation_type, "joint")) {
      obj <- object
      obj$truncation <- NULL
      j <- jt[1]; Fa <- rep_len(Fab[[1]]$Fa, nrow(U)); Fb <- rep_len(Fab[[1]]$Fb, nrow(U))
      core <- function(uj) {
        Uj <- U
        Uj[, j] <- .tc_clip(uj)
        .tc_jointprob_u(obj, Uj, newdata)
      }
      hi <- core(pmin(pmax(U[, j], Fa), Fb))
      lo <- ifelse(Fa > 1e-12, core(Fa), 0)
      return(pmin(pmax((hi - lo) / (Fb - Fa), 0), 1))
    }
    for (k in seq_along(jt)) {
      j <- jt[k]
      U[, j] <- .tc_clip(pmin(pmax((U[, j] - Fab[[k]]$Fa) / (Fab[[k]]$Fb - Fab[[k]]$Fa), 0), 1))
    }
  }
  .tc_jointprob_u(object, U, newdata)
}

# Copula CDF at the rows of U (n x 3 uniforms) for the dependence parameters
# of each row of newdata.
.tc_jointprob_u <- function(object, U, newdata = NULL) {
  pcc <- object$copula == "pcc"
  U <- cbind(matrix(U, ncol = 3),
             if (pcc) do.call(cbind, .tc_new_dep_eta(object, newdata)) else
               .tc_new_rho(object, newdata)[, 1:3, drop = FALSE])
  gl <- .tc_gauss_legendre(48L)
  ptau <- if (pcc) .tc_pcc_tau(object, newdata)
  # covariate-dependent df: one more column per df predictor (columns 7 on)
  vdf <- isFALSE(object$df_constant) && length(object$index$df) > 0L
  if (vdf) {
    tl <- if (pcc) ptau else .tc_new_tau(object, newdata)
    U <- cbind(U, matrix(unlist(lapply(tl, function(t) {
      if (is.null(t)) rep(NA_real_, nrow(U)) else rep_len(t, nrow(U))
    })), nrow = nrow(U)))
  }
  ord <- if (pcc) .tc_pcc_order(object$root)
  key <- apply(U, 1, paste, collapse = "_")
  first <- !duplicated(key)
  val <- apply(U[first, , drop = FALSE], 1, function(r) {
    if (pcc) {
      pt <- if (!vdf) ptau else lapply(1:3, function(k) if (is.na(r[[6 + k]])) NULL else r[[6 + k]])
      return(.tc_pcc_cdf(r[ord], r[4:6], object$families, pt, gl, object$nonsimplified,
                         object$coefficients[object$index$ns]))
    }
    .tc_copula_cdf(r[1:3], .tc_R(r[[4]], r[[5]], r[[6]]),
                   if (vdf) 2.05 + exp(r[[7]]) else object$df, object$copula, gl)
  })
  unname(val[match(key, key[first])])
}

#' Estimation and specification diagnostics
#'
#' Reports optimizer convergence, the size of the gradient at the optimum,
#' eigenvalues and condition number of the observed information, the
#' log-likelihood gain of full ML over the two-stage (IFM) start, the
#' agreement of all optimizer starts, Kolmogorov-Smirnov tests of uniformity
#' of each marginal PIT, and the largest gap between the fitted copula
#' correlation matrix and the empirical correlation of the copula scores.
#' The KS p-values ignore parameter estimation and are therefore
#' conservative screening tools, not exact tests.
#'
#' @param object A `tricopula` fit.
#' @return An object of class `tc_diagnostics` (a list).
#' @export
tc_diagnostics <- function(object) {
  .tc_check(object)
  ev <- if (!is.null(object$hessian)) {
    eigen(object$hessian, symmetric = TRUE, only.values = TRUE)$values
  } else NA_real_
  U <- object$uniforms
  # PCC: the Rosenblatt transforms (u1, h12, h23|1) are independent uniforms
  # under a correct model; report the largest correlation of their normal
  # scores (target 0) in place of the score-correlation gap.
  Z <- if (object$copula == "pcc") NULL else
    if (object$copula == "gaussian") stats::qnorm(U) else
      stats::qt(U, if (is.null(object$df_i)) object$df else object$df_i)
  gap <- if (object$copula == "pcc") {
    C <- stats::cor(stats::qnorm(object$rosenblatt), use = "complete.obs")
    max(abs(C[upper.tri(C)]))
  } else max(abs(stats::cor(Z, use = "pairwise.complete.obs") - object$correlation))
  # For the t copula the scores' correlation equals R only if df > 2.
  ks <- vapply(1:3, function(j) {
    suppressWarnings(stats::ks.test(U[, j], "punif")$p.value)
  }, numeric(1))
  names(ks) <- object$response_names
  structure(list(
    convergence = object$convergence,
    newton_iterations = object$newton_iterations,
    max_abs_gradient = max(abs(object$gradient)),
    hessian_min_eigenvalue = min(ev), hessian_max_eigenvalue = max(ev),
    hessian_condition_number = max(ev) / min(ev),
    hessian_positive_definite = isTRUE(min(ev) > 0),
    logLik = object$logLik,
    logLik_ifm = if (is.null(object$ifm)) NA_real_ else object$ifm$logLik,
    start_logLik_range = diff(range(object$starts$logLik)),
    pit_ks_pvalue = ks,
    max_abs_score_correlation_gap = gap
  ), class = "tc_diagnostics")
}

#' @export
print.tc_diagnostics <- function(x, ...) {
  cat("TriCopula diagnostics\n")
  cat("  optimizer convergence code     :", x$convergence, "\n")
  cat("  Newton refinement steps        :", x$newton_iterations, "\n")
  cat("  max |gradient| at optimum      :", format(x$max_abs_gradient, digits = 3), "\n")
  cat("  information positive definite  :", x$hessian_positive_definite, "\n")
  cat("  information condition number   :", format(x$hessian_condition_number, digits = 4), "\n")
  cat("  logLik (full ML)               :", format(x$logLik, digits = 8), "\n")
  cat("  logLik (two-stage IFM)         :", format(x$logLik_ifm, digits = 8), "\n")
  cat("  logLik range across starts     :", format(x$start_logLik_range, digits = 3), "\n")
  cat("  PIT uniformity KS p-values     :", paste(names(x$pit_ks_pvalue),
      format(x$pit_ks_pvalue, digits = 3), sep = " = ", collapse = ", "), "\n")
  cat("  max |cor(scores) - R|          :", format(x$max_abs_score_correlation_gap, digits = 3), "\n")
  invisible(x)
}

#' Likelihood-ratio test of two nested fits
#'
#' Typical use: constant dependence (`dependence = ~ 1`) against
#' covariate-dependent dependence. Nesting is the user's responsibility. The
#' Gaussian copula is the `df -> Inf` boundary of the Student-t copula, so the
#' chi-square reference distribution is conservative for that comparison.
#'
#' @param restricted,full Two `tricopula` fits to the same data.
#' @return A data frame with the statistic, degrees of freedom and p-value.
#' @export
tc_lrtest <- function(restricted, full) {
  .tc_check(restricted); .tc_check(full)
  if (restricted$n != full$n) stop("The fits use different numbers of observations.")
  if (!is.null(restricted$penalties) || !is.null(full$penalties)) {
    warning("At least one fit has penalized s() terms: the chi-square reference with ",
            "the difference in parameter counts is not valid; compare by AIC.",
            call. = FALSE)
  }
  k <- full$parameter_count - restricted$parameter_count
  if (k <= 0) stop("'full' must have more parameters than 'restricted'.")
  stat <- 2 * (full$logLik - restricted$logLik)
  data.frame(logLik_restricted = restricted$logLik, logLik_full = full$logLik,
             statistic = stat, df = k,
             p.value = stats::pchisq(stat, k, lower.tail = FALSE))
}

#' Compare fitted models by information criteria
#'
#' @param ... Two or more `tricopula` fits to the same responses.
#' @return A data frame with log-likelihood, parameter count, AIC and BIC.
#' @export
tc_compare <- function(...) {
  fits <- list(...)
  lapply(fits, .tc_check)
  nm <- vapply(as.list(substitute(list(...)))[-1L], deparse, character(1))
  data.frame(model = nm,
             copula = vapply(fits, function(f) {
               if (f$copula == "pcc") paste0("pcc(", paste(f$families, collapse = ","), ")") else
                 f$copula
             }, character(1)),
             margins = vapply(fits, function(f) paste(f$margins, collapse = "/"), character(1)),
             logLik = vapply(fits, `[[`, numeric(1), "logLik"),
             npar = vapply(fits, `[[`, integer(1), "parameter_count"),
             AIC = vapply(fits, stats::AIC, numeric(1)),
             BIC = vapply(fits, stats::BIC, numeric(1)))
}
