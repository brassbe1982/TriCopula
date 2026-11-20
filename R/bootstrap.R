# Bootstrap inference ------------------------------------------------------

# Refit the model of `object` on `data`, started from its estimates. Smooth
# terms reuse the spline bases of `object`; `sp` gives the smoothing
# parameters (default: those of `object`, held fixed; NULL selects them
# again). Returns the fit or a try-error.
.tc_refit <- function(object, data, sp = object$sp, hessian = FALSE,
                      start = object$coefficients, weights = object$weights,
                      simplified = !isFALSE(object$simplified),
                      k_cond = if (is.null(object$nonsimplified)) 6L else
                        object$nonsimplified$k) {
  try(tricopula(object$formulas, data, margins = object$margins, sp = sp,
                weights = weights,
                weights_type = if (is.null(object$weights_type)) "auto" else object$weights_type,
                missing = if (identical(object$missing, "likelihood")) "likelihood" else "fail",
                copula = object$copula, start = start,
                control = object$control, hessian = hessian,
                dependence = object$dependence,
                scale = object$scale_formulas,
                families = object$families,
                df = if (isFALSE(object$df_constant)) {
                  if (object$copula == "student") object$df_formulas[[1]] else
                    object$df_formulas
                } else ~ 1,
                root = if (is.null(object$root)) 1L else object$root,
                shape = if (is.null(object$shape_design)) ~ 1 else
                  object$shape_formulas,
                smooth_basis = if (!is.null(object$smooth_blocks)) object$smooths,
                simplified = simplified, k_cond = k_cond,
                truncation = object$truncation,
                truncation_type = if (is.null(object$truncation_type)) "margin" else
                  object$truncation_type,
                censoring = object$censoring,
                cens_nodes = if (is.null(object$cens_nodes)) c(64L, 32L, 16L) else
                  object$cens_nodes,
                cens_method = if (is.null(object$cens_method)) "quadrature" else
                  object$cens_method,
                cens_check = FALSE),
      silent = TRUE)
}

.tc_boot_refit <- function(d, object, reselect = FALSE) {
  w <- attr(d, "tc_weights")
  if (is.null(w)) w <- object$weights
  f <- if (reselect) .tc_refit(object, d, sp = NULL, weights = w) else
    .tc_refit(object, d, weights = w)
  k <- length(object$coefficients) + length(.tc_natural(object$coefficients, object))
  if (inherits(f, "try-error") || f$convergence != 0) return(rep(NA_real_, k))
  c(f$coefficients, .tc_natural(f$coefficients, f))
}

# Quantities that are defined whatever copula type, pair-copula families,
# margins or covariates were selected (0.14.1): Kendall's tau and the
# tail-dependence coefficients of the pairs, the average conditional means
# over the rows of `newdata`, and average joint probabilities
# P(Y1 <= t1, Y2 <= t2, Y3 <= t3 | x) for each row of `thresholds`.
.tc_common <- function(par, object, thresholds = NULL, newdata = NULL,
                       what = c("dep", "mean", "jointprob")) {
  object$coefficients[] <- par
  out <- numeric(0)
  if ("dep" %in% what) {
    nat <- .tc_natural(par, object)
    if (object$copula == "pcc") {
      pl <- .tc_pcc_labels(object$root)$short
      td <- tc_taildep(object)
      out <- c(nat[paste0("tau", pl)],
               stats::setNames(c(td$lower, td$upper),
                               c(paste0("lambdaL", pl[1:2]), paste0("lambdaU", pl[1:2]))))
    } else {
      r <- unname(nat[c("rho12", "rho13", "rho23")])
      pc <- function(a, b, c) (a - b * c) / sqrt((1 - b^2) * (1 - c^2))
      rr <- c(r, unname(nat[["rho23|1"]]), pc(r[2], r[1], r[3]), pc(r[1], r[2], r[3]))
      lam <- if (object$copula == "gaussian") rep(0, 3) else {
        nu <- unname(nat[["df"]])
        2 * stats::pt(-sqrt((nu + 1) * (1 - r) / (1 + r)), nu + 1)
      }
      out <- stats::setNames(
        c(2 / pi * asin(rr), lam, lam),
        c("tau12", "tau13", "tau23", "tau23|1", "tau13|2", "tau12|3",
          paste0("lambdaL", c("12", "13", "23")), paste0("lambdaU", c("12", "13", "23"))))
    }
  }
  if ("mean" %in% what) {
    out <- c(out, stats::setNames(colMeans(.tc_quantity(par, object, newdata, "mean", NULL)),
                                  paste0("mean", 1:3)))
  }
  if ("jointprob" %in% what && !is.null(thresholds)) {
    jp <- vapply(seq_len(nrow(thresholds)), function(i) {
      mean(.tc_quantity(par, object, newdata, "jointprob", thresholds[i, ]))
    }, numeric(1))
    out <- c(out, stats::setNames(jp, paste0("jointprob", seq_len(nrow(thresholds)))))
  }
  out
}

# Estimates and delta-method standard errors of .tc_common() for a fit.
.tc_common_table <- function(object, thresholds = NULL, newdata = NULL) {
  b <- object$coefficients
  est <- .tc_common(b, object, thresholds, newdata)
  se <- stats::setNames(rep(NA_real_, length(est)), names(est))
  if (!is.null(object$vcov)) {
    V <- stats::vcov(object)
    V[is.na(V)] <- 0
    sef <- function(J) sqrt(pmax(diag(J %*% V %*% t(J)), 0))
    J1 <- numDeriv::jacobian(.tc_common, b, object = object, newdata = newdata,
                             what = c("dep", "mean"))
    se[seq_len(nrow(J1))] <- sef(J1)
    if (!is.null(thresholds)) {
      J2 <- numDeriv::jacobian(.tc_common, b, method = "simple",
                               method.args = list(eps = 1e-4), object = object,
                               thresholds = thresholds, newdata = newdata,
                               what = "jointprob")
      se[nrow(J1) + seq_len(nrow(J2))] <- sef(J2)
    }
  }
  data.frame(estimate = unname(est), std.error = unname(se), row.names = names(est))
}

# One bootstrap replicate with the whole selection repeated (0.14.0; margins,
# covariates and copula type since 0.14.1). Returns NULL if the selection or
# the fit fails.
.tc_boot_select <- function(d, object, thresholds = NULL, newdata = NULL) {
  s <- object$selection
  w <- attr(d, "tc_weights")
  if (is.null(w)) w <- object$weights
  dots <- s$dots
  dots$hessian <- FALSE
  dots$start <- NULL
  or <- function(a, b) if (is.null(a)) b else a
  r <- try(suppressWarnings(do.call(tc_select, c(list(
    or(s$formulas, object$formulas), d, margins = or(s$margins, object$margins),
    candidates = s$candidates, criterion = s$criterion,
    scale = or(s$scale, object$scale_formulas),
    dependence = or(s$dependence, object$dependence), root = s$root,
    shape = or(s$shape, if (is.null(object$shape_design)) ~ 1 else object$shape_formulas),
    weights = w,
    missing = if (identical(object$missing, "likelihood")) "likelihood" else "fail",
    copula = or(s$copula, "pcc"), covariates = or(s$covariates, "none")),
    dots))), silent = TRUE)
  if (inherits(r, "try-error") || r$fit$convergence != 0) return(NULL)
  f <- r$fit
  cm <- try(suppressWarnings(.tc_common(f$coefficients, f, thresholds, newdata)), silent = TRUE)
  if (inherits(cm, "try-error")) return(NULL)
  list(families = r$families, root = if (is.null(r$root)) NA_integer_ else r$root,
       copula = f$copula, margins = f$margins,
       terms = lapply(f$formulas, function(g) attr(stats::terms(g), "term.labels")),
       values = c(f$coefficients, .tc_natural(f$coefficients, f)), common = cm)
}

# Bootstrap with selection: summary of the replicates. The main tables hold
# quantities that mean the same under every selection: the coefficients and
# natural parameters of the margins whose family was not selected, Kendall's
# tau of the pairs and (`derived`) tail dependence, average means and joint
# probabilities. `conditional` holds the copula parameters, from the
# replicates that made the same choice as the original fit.
.tc_boot_select_summary <- function(object, res, B, type, level, common0 = NULL) {
  ok <- !vapply(res, is.null, logical(1))
  if (sum(ok) < 2L) stop("Fewer than two bootstrap replicates could be selected and fitted.")
  res <- res[ok]
  R <- length(res)
  s <- object$selection
  nat <- .tc_natural(object$coefficients, object)
  est <- c(object$coefficients, nat)
  msel <- if (is.list(s$margins)) lengths(s$margins) > 1L else rep(FALSE, 3L)
  cov_sel <- identical(s$covariates, "backward")
  nm_m <- names(est)[grepl("^margin", names(est))]
  if (cov_sel) {
    # every coefficient of the full mean formulas is reported
    full <- unlist(lapply(1:3, function(j) {
      tt <- stats::delete.response(stats::terms(s$formulas[[j]]))
      paste0("margin", j, "_", colnames(stats::model.matrix(tt, object$data)))
    }))
    nm_m <- union(union(nm_m, full), unlist(lapply(res, function(r) {
      grep("^margin", names(r$values), value = TRUE)
    })))
    nm_m <- nm_m[order(as.integer(substr(nm_m, 7L, 7L)), nm_m %in% names(nat))]
  }
  if (any(msel)) {
    nm_m <- nm_m[!grepl(paste0("^margin(", paste(which(msel), collapse = "|"), ")_"), nm_m)]
  }
  keep <- c(nm_m, names(est)[grepl("^tau", names(est))])
  pick <- function(v) {
    x <- stats::setNames(unname(v[keep]), keep)
    # a term that was not selected has coefficient zero
    if (cov_sel) x[is.na(x) & grepl("^margin", keep)] <- 0
    x
  }
  M <- do.call(rbind, lapply(res, function(r) unname(pick(r$values))))
  colnames(M) <- keep
  a <- (1 - level) / 2
  qf <- function(x, p) if (sum(!is.na(x)) < 2L) NA_real_ else
    unname(stats::quantile(x, probs = p, na.rm = TRUE))
  sdf <- function(x) if (sum(!is.na(x)) < 2L) NA_real_ else stats::sd(x, na.rm = TRUE)
  mk <- function(e, Mx) {
    data.frame(estimate = unname(e), boot.se = apply(Mx, 2, sdf),
               lower = apply(Mx, 2, qf, p = a), upper = apply(Mx, 2, qf, p = 1 - a),
               replicates = colSums(!is.na(Mx)), row.names = colnames(Mx))
  }
  tab <- mk(pick(est), M)
  tab$model.se <- NA_real_
  dn <- NULL
  if (!is.null(object$vcov)) {
    dn <- tc_dependence(object)
    se <- c(sqrt(diag(stats::vcov(object))), stats::setNames(dn$std.error, rownames(dn)))
    tab$model.se <- unname(se[keep])
  }
  # derived quantities (0.14.1)
  derived <- NULL; MD <- NULL
  if (!is.null(common0)) {
    dnm <- setdiff(rownames(common0), keep)
    MD <- do.call(rbind, lapply(res, function(r) unname(r$common[dnm])))
    colnames(MD) <- dnm
    derived <- mk(common0[dnm, "estimate"], MD)
    derived$model.se <- common0[dnm, "std.error"]
  }
  # copula parameters, conditional on the same choice (0.14.1)
  pcc0 <- object$copula == "pcc"
  root0 <- if (is.null(object$root)) 1L else object$root
  lab0 <- .tc_pcc_labels(root0)$short
  cn <- names(nat)[grepl("^(par|df|rho)", names(nat))]
  same_choice <- function(r, nm) {
    if (!pcc0) return(identical(r$copula, object$copula))
    k <- match(sub("^(par|df)", "", nm), lab0)
    identical(r$copula, "pcc") && isTRUE(r$root == root0) &&
      identical(r$families[k], object$families[k])
  }
  MC <- vapply(cn, function(nm) {
    vapply(res, function(r) if (same_choice(r, nm)) unname(r$values[nm]) else NA_real_,
           numeric(1))
  }, numeric(R))
  MC <- matrix(MC, R, length(cn), dimnames = list(NULL, cn))
  conditional <- mk(nat[cn], MC)
  conditional$model.se <- if (is.null(dn)) NA_real_ else dn[cn, "std.error"]
  # what was selected
  is_pcc <- vapply(res, function(r) identical(r$copula, "pcc"), logical(1))
  roots <- vapply(res, function(r) as.integer(r$root), integer(1))
  sel <- data.frame(pair = character(0), family = character(0), share = numeric(0),
                    original = logical(0))
  if (any(is_pcc)) {
    fams <- do.call(rbind, lapply(res[is_pcc], `[[`, "families"))
    labs <- t(vapply(roots[is_pcc], function(r) .tc_pcc_labels(r)$short, character(3)))
    sel <- as.data.frame(table(pair = as.vector(labs), family = as.vector(fams)),
                         stringsAsFactors = FALSE)
    sel <- sel[sel$Freq > 0, ]
    sel$share <- sel$Freq / R
    sel$Freq <- NULL
    sel$original <- pcc0 & paste(sel$pair, sel$family) %in% paste(lab0, object$families)
    sel <- sel[order(match(sel$pair, lab0), -sel$share), ]
    rownames(sel) <- NULL
  }
  share_tab <- function(x, original, what) {
    t0 <- table(x) / R
    d <- data.frame(names(t0), as.numeric(t0), names(t0) == original,
                    stringsAsFactors = FALSE)
    names(d) <- c(what, "share", "original")
    d[order(-d$share), , drop = FALSE]
  }
  copulas <- if (length(s$copula) > 1L) {
    share_tab(vapply(res, `[[`, "", "copula"), object$copula, "copula")
  }
  margin_selection <- if (any(msel)) {
    do.call(rbind, lapply(which(msel), function(j) {
      cbind(margin = j, share_tab(vapply(res, function(r) r$margins[j], ""),
                                  object$margins[j], "family"))
    }))
  }
  terms0 <- lapply(object$formulas, function(g) attr(stats::terms(g), "term.labels"))
  covariate_selection <- if (cov_sel) {
    do.call(rbind, lapply(1:3, function(j) {
      full <- attr(stats::terms(s$formulas[[j]]), "term.labels")
      if (!length(full)) return(NULL)
      data.frame(margin = j, term = full,
                 share = vapply(full, function(t) {
                   mean(vapply(res, function(r) t %in% r$terms[[j]], logical(1)))
                 }, numeric(1)),
                 original = full %in% terms0[[j]], row.names = NULL)
    }))
  }
  same <- mean(vapply(res, function(r) {
    identical(r$copula, object$copula) && identical(r$margins, object$margins) &&
      identical(r$terms, terms0) &&
      (!pcc0 || (isTRUE(r$root == root0) && identical(r$families, object$families)))
  }, logical(1)))
  is_nat <- keep %in% names(nat)
  structure(list(coefficients = tab[!is_nat, ], natural = tab[is_nat, ],
                 replicates = M, failures = sum(!ok), B = B, type = type, level = level,
                 select = TRUE, selection = sel,
                 roots = table(root = factor(roots, levels = 1:3)) / R,
                 same_selection = same, derived = derived, derived_replicates = MD,
                 conditional = conditional, copulas = copulas,
                 margin_selection = margin_selection,
                 covariate_selection = covariate_selection), class = "tc_bootstrap")
}

#' Bootstrap standard errors and confidence intervals
#'
#' Refits the model on bootstrap samples, starting each fit from the original
#' estimates. The parametric bootstrap simulates new responses from the
#' fitted model (covariates held fixed); the nonparametric bootstrap resamples
#' rows with replacement and remains valid when the model is misspecified.
#'
#' @param object A `tricopula` fit.
#' @param B Number of bootstrap samples.
#' @param type `"parametric"` or `"nonparametric"`.
#' @param level Confidence level of the percentile intervals.
#' @param seed Optional random seed (samples are drawn in the calling
#'   session, so results do not depend on `cl`).
#' @param cl Optional cluster from [parallel::makeCluster()] on which the
#'   refits are run; the package must be installed on the workers.
#' @return An object of class `tc_bootstrap`: a list with `coefficients` and
#'   `natural` (data frames with the estimate, bootstrap standard error,
#'   percentile interval and, for comparison, the model-based standard
#'   error), the matrix of `replicates`, and the number of `failures`
#'   (non-converged refits, which are dropped).
#' @param sp For fits with penalized smooth terms: `"fixed"` (the default)
#'   holds the smoothing parameters at their estimates in every refit;
#'   `"select"` selects them again in every refit, which is slower but lets
#'   the bootstrap reflect their uncertainty.
#' @param select Repeat the selection of the pair-copula families (and of
#'   the root, if it was selected) in every bootstrap sample? Requires the
#'   fit returned by [tc_select()]. The intervals then include the
#'   uncertainty of the selection. Because the dependence parameters mean
#'   different things under different families, only quantities common to
#'   all families are reported: the coefficients and natural parameters of
#'   the margins and Kendall's tau of the pairs. The result also holds
#'   `selection` (how often each family was selected for each pair),
#'   `roots`, and `same_selection` (the share of samples in which the
#'   selection of the original fit was reproduced). A tau whose pair does
#'   not exist in a replicate with another root is missing there; column
#'   `replicates` counts the values used.
#'
#'   Since 0.14.1 every step of [tc_select()] is repeated (margins,
#'   covariates, copula type, families, root) and the result also holds
#'   `derived`: Kendall's tau where not already in `natural`, the lower and
#'   upper tail-dependence coefficients of the unconditional pairs
#'   (`lambdaL12`, `lambdaU12`, ...), the conditional means averaged over
#'   the rows of `newdata` (`mean1` to `mean3`) and, if `thresholds` is
#'   given, joint probabilities (`jointprob1`, ...); `conditional`: the
#'   copula parameters, computed from those samples only in which the same
#'   family (or copula type) was chosen as in the original fit, which are
#'   intervals given that choice and not intervals after selection; and
#'   `copulas`, `margin_selection`, `covariate_selection` (selection
#'   frequencies), where these were selected. Coefficients of a margin whose
#'   family was selected are not reported (they mean different things under
#'   different families; use the means). With `covariates = "backward"` a
#'   term that was not selected counts as a coefficient of zero.
#' @param thresholds With `select = TRUE`: three thresholds `(t1, t2, t3)`,
#'   or a matrix with three columns, for joint probabilities
#'   `P(Y1 <= t1, Y2 <= t2, Y3 <= t3 | x)` averaged over the rows of
#'   `newdata` (see [tc_jointprob()]).
#' @param newdata With `select = TRUE`: covariate rows at which the means
#'   and joint probabilities are evaluated; defaults to the estimation data.
#'   Joint probabilities cost one trivariate integral per distinct row and
#'   bootstrap sample, so pass a small `newdata` with `thresholds`.
#' @details Prior weights are kept with their rows. With missing responses
#'   (`missing = "likelihood"`) the parametric bootstrap simulates complete
#'   responses and then deletes those that are missing in the data, which
#'   assumes that the missingness does not depend on the responses.
#'
#'   The parametric bootstrap requires each response to be a plain
#'   column of the data (not a transformation inside the formula). For fits
#'   with smooth terms both bootstraps keep the spline bases (knots) of the
#'   original fit, so the coefficients of the replicates are comparable.
#' @export
tc_bootstrap <- function(object, B = 200L, type = c("parametric", "nonparametric"),
                         level = 0.95, seed = NULL, cl = NULL,
                         sp = c("fixed", "select"), select = FALSE,
                         thresholds = NULL, newdata = NULL) {
  .tc_check(object)
  type <- match.arg(type)
  sp <- match.arg(sp)
  if (is.null(object$data)) stop("The fit does not contain the estimation data.")
  if (!is.logical(select) || length(select) != 1L || is.na(select)) {
    stop("select must be TRUE or FALSE.")
  }
  if (select) {
    if (is.null(object$selection)) {
      stop("select = TRUE needs the fit returned by tc_select() (its element `fit`).")
    }
    if (!is.null(object$selection$split)) {
      stop("select = TRUE is not available for a fit from tc_select(split = ): ",
           "its standard errors are already valid given the selection.")
    }
    if (!is.null(object$smooth_blocks)) {
      stop("select = TRUE is not available for fits with smooth terms.")
    }
  }
  common0 <- NULL
  if (!select && (!is.null(thresholds) || !is.null(newdata))) {
    stop("thresholds and newdata are only used with select = TRUE.")
  }
  if (select) {
    if (!is.null(thresholds)) {
      if (!is.numeric(thresholds) || (is.null(dim(thresholds)) && length(thresholds) != 3L) ||
          (!is.null(dim(thresholds)) && ncol(thresholds) != 3L) || anyNA(thresholds)) {
        stop("thresholds must be three numbers or a matrix with three columns.")
      }
      thresholds <- matrix(thresholds, ncol = 3L)
    }
    if (is.null(newdata)) newdata <- object$data
    common0 <- .tc_common_table(object, thresholds, newdata)
  }
  if (type == "parametric" && !is.null(object$censored)) {
    stop("The parametric bootstrap is not available for fits with censored responses ",
         "(the censoring mechanism is not modelled); use type = \"nonparametric\".")
  }
  if (select && (!is.null(object$censored) || !is.null(object$truncation))) {
    stop("select = TRUE is not available for fits with truncation or censoring.")
  }
  if (!is.null(seed)) set.seed(seed)
  rn <- object$response_names
  if (type == "parametric" && !all(rn %in% names(object$data))) {
    stop("The parametric bootstrap needs responses that are plain columns of the data.")
  }
  n <- nrow(object$data)
  samples <- lapply(seq_len(B), function(b) {
    d <- object$data
    if (type == "parametric") {
      sim <- stats::simulate(object)
      # responses that are missing in the data stay missing
      if (!is.null(object$observed)) sim[!object$observed] <- NA
      d[rn] <- sim
    } else {
      ii <- sample.int(n, replace = TRUE)
      d <- d[ii, , drop = FALSE]
      if (!is.null(object$weights)) attr(d, "tc_weights") <- object$weights[ii]
    }
    d
  })
  if (select) {
    res <- if (is.null(cl)) {
      lapply(samples, .tc_boot_select, object = object, thresholds = thresholds,
             newdata = newdata)
    } else {
      parallel::parLapply(cl, samples, .tc_boot_select, object = object,
                          thresholds = thresholds, newdata = newdata)
    }
    return(.tc_boot_select_summary(object, res, B, type, level, common0))
  }
  res <- if (is.null(cl)) {
    lapply(samples, .tc_boot_refit, object = object, reselect = sp == "select")
  } else {
    parallel::parLapply(cl, samples, .tc_boot_refit, object = object,
                        reselect = sp == "select")
  }
  M <- do.call(rbind, res)
  nat <- .tc_natural(object$coefficients, object)
  colnames(M) <- c(names(object$coefficients), names(nat))
  ok <- stats::complete.cases(M)
  if (sum(ok) < 2L) stop("Fewer than two bootstrap refits converged.")
  M <- M[ok, , drop = FALSE]
  a <- (1 - level) / 2
  tab <- data.frame(estimate = c(object$coefficients, nat),
                    boot.se = apply(M, 2, stats::sd),
                    lower = apply(M, 2, stats::quantile, probs = a),
                    upper = apply(M, 2, stats::quantile, probs = 1 - a),
                    row.names = colnames(M))
  p <- length(object$coefficients)
  tab$model.se <- NA_real_
  if (!is.null(object$vcov)) {
    tab$model.se[seq_len(p)] <- sqrt(diag(stats::vcov(object)))
    tab$model.se[-seq_len(p)] <- tc_dependence(object)$std.error
  }
  structure(list(coefficients = tab[seq_len(p), ], natural = tab[-seq_len(p), ],
                 replicates = M, failures = sum(!ok), B = B, type = type,
                 level = level), class = "tc_bootstrap")
}

#' @export
print.tc_bootstrap <- function(x, ...) {
  cat("TriCopula", x$type, "bootstrap:", x$B, "samples,", x$failures,
      "failed refits,", 100 * x$level, "% percentile intervals\n\n")
  print(round(x$coefficients, 4))
  cat("\nNatural scale:\n")
  print(round(x$natural, 4))
  if (isTRUE(x$select)) {
    cat("\nSelection repeated in every sample; the selection of the original fit was\n",
        "reproduced in ", format(100 * x$same_selection, digits = 3),
        "% of them. Selected families by pair:\n", sep = "")
    s <- x$selection
    s$share <- round(s$share, 3)
    print(s, row.names = FALSE)
    show <- function(title, d) {
      if (is.null(d) || !nrow(d)) return(invisible())
      cat("\n", title, "\n", sep = "")
      num <- vapply(d, is.numeric, logical(1))
      d[num] <- lapply(d[num], round, 4)
      print(d, row.names = !is.null(d$estimate))
    }
    show("Copula type (share of the samples):", x$copulas)
    show("Margins (share of the samples):", x$margin_selection)
    show("Terms of the mean formulas (share of the samples that kept them):",
         x$covariate_selection)
    show("Tail dependence, average means and joint probabilities:", x$derived)
    show(paste("Copula parameters, from the samples with the same choice as the",
               "original fit\n(conditional on that choice; see 'replicates'):"),
         x$conditional)
  }
  invisible(x)
}
