# Pair-copula family selection ----------------------------------------------

# Maximum-likelihood fit of one pair copula to pseudo-observations (u, v).
.tc_fit_pair <- function(u, v, family) {
  is_t <- .tc_pair_base(family) == "student"
  kt <- stats::cor(u, v, method = "kendall")
  st <- c(.tc_pair_eta_from_tau(kt, family), if (is_t) log(8 - 2.05))
  fn <- function(th) {
    l <- try(suppressWarnings(.tc_pair(u, v, th[1], family, if (is_t) th[2])$l),
             silent = TRUE)
    if (inherits(l, "try-error")) return(1e100)
    s <- -sum(l)
    if (is.finite(s)) s else 1e100
  }
  opt <- if (is_t) {
    stats::optim(st, fn, method = "BFGS", control = list(maxit = 300))
  } else {
    lim <- .tc_pair_eta_lim(family)
    o <- stats::optimize(fn, c(-lim, lim))
    o2 <- stats::optim(st, fn, method = "BFGS", control = list(maxit = 300))
    if (o$objective < o2$value) list(par = o$minimum, value = o$objective) else o2
  }
  list(eta = opt$par[1], tau = if (is_t) opt$par[2], logLik = -opt$value,
       npar = length(st))
}

#' Select the pair-copula families of a PCC and fit the selected model
#'
#' Chooses the three families of the C-vine `c12`, `c13`, `c23|1` and then
#' fits the selected model by joint maximum likelihood.
#'
#' @details
#' Selection is sequential, as is usual for vines: the margins are fitted by
#' marginal maximum likelihood; each candidate family is fitted to the
#' pseudo-observations of pair (1, 2) and of pair (1, 3) and the best one by
#' `criterion` is kept; the candidates are then fitted to the h-function
#' transforms for the conditional pair (2, 3 | 1). Dependence is taken as
#' constant during selection. Only the selected combination is fitted
#' jointly, with `dependence` and `scale` as supplied.
#'
#' By default the candidates for a pair with positive Kendall's tau are
#' `gaussian`, `student`, `clayton`, `gumbel`, `clayton180`, `gumbel180`;
#' for negative tau the Clayton and Gumbel families are replaced by their
#' 90- and 270-degree rotations. `candidates = "all"` adds `frank`, `joe`
#' and `joe180` (for negative tau: `frank`, `joe90`, `joe270`).
#'
#' With `root = 2` or `root = 3` the pairs are those of the vine rooted at
#' that response (see [tricopula()]). With `root = "auto"` the sequential
#' selection and the joint fit are carried out for each of the three roots
#' and the root whose joint fit has the smallest `criterion` is returned;
#' this takes three joint fits.
#'
#' The standard errors of the final fit do not account for the selection of
#' the families or of the root. Two remedies (0.14.0): with `split` the
#' families (and the root) are selected on a random part of the rows and the
#' selected model is estimated on the remaining rows only, so that its
#' standard errors are valid given the selection, at the price of a smaller
#' estimation sample; `tc_bootstrap(fit, select = TRUE)` repeats the
#' selection in every bootstrap sample.
#'
#' **Selection of the margins, the covariates and the copula type (0.14.1).**
#' If `margins` is a list of three character vectors, the family of each
#' margin is chosen first, among its candidates, by the criterion of the
#' marginal maximum-likelihood fit. With `covariates = "backward"` the terms
#' of each mean formula are then removed one at a time, as long as that
#' lowers the criterion of the marginal fit (terms of the `scale`, `shape`
#' and `dependence` formulas are not selected). The pair-copula families are
#' selected with the chosen margins and formulas. If `copula` names more
#' than one copula type, each type is fitted jointly and the one with the
#' smallest criterion of the joint fit is returned. Sample splitting and
#' `tc_bootstrap(select = TRUE)` repeat all of these steps.
#'
#' @param formulas,data,scale,dependence,shape,weights,missing As in
#'   [tricopula()]. The selection of the pair-copula families uses the
#'   complete pairs and no weights; the marginal selection steps and the
#'   final joint fit use both.
#' @param margins Three margin names as in [tricopula()], or a list of three
#'   character vectors of candidate names, e.g.
#'   `list("gaussian", c("lognormal", "gamma"), c("gamma", "weibull"))`.
#' @param copula The copula types considered: `"pcc"` (the default, as in
#'   earlier versions), any of `"gaussian"`, `"student"`, `"pcc"`, or
#'   `"auto"` for all three.
#' @param covariates `"none"` (the default) or `"backward"`: backward
#'   elimination of the terms of the three mean formulas.
#' @param candidates Optional character vector of candidate families
#'   (applied to all three pairs), or `"all"`.
#' @param criterion `"AIC"` or `"BIC"`.
#' @param root Root of the vine: 1 (the default), 2, 3 or `"auto"`.
#' @param split `NULL` (the default: selection and estimation on all rows) or
#'   a number strictly between 0 and 1: the share of the rows, drawn at
#'   random, that is used for the selection only. The returned `fit` is then
#'   estimated on the other rows.
#' @param seed Optional random seed for the split.
#' @param ... Further arguments passed to [tricopula()] for the final fit.
#' @return A list with `families` (the selection), `table` (log-likelihood
#'   and criterion of every candidate for every pair), `fit` (the joint
#'   `tricopula` fit of the selected PCC) and `root`. With `root = "auto"`
#'   also `roots`, a data frame with the selected families, log-likelihood
#'   and criterion of the joint fit for each root; `table` then refers to
#'   the selected root. With `split` also `rows_select` and `rows_estimate`
#'   (row numbers of `data`) and `selection_fit` (the joint fit on the
#'   selection part, without covariance matrix); `table` refers to the
#'   selection part. `fit$selection` holds the settings of the selection.
#'   Since 0.14.1 also `copula`, `margins` and `formulas` (those of the
#'   returned fit) and, if these were selected, `copulas` (log-likelihood
#'   and criterion of the joint fit per copula type), `margin_table` and
#'   `covariate_table`. If an elliptical copula is returned, `families` and
#'   `root` are `NULL` and `pcc_families` holds the families of the PCC that
#'   was compared.
#' @export
tc_select <- function(formulas, data,
                      margins = c("gaussian", "lognormal", "gamma"),
                      candidates = NULL, criterion = c("AIC", "BIC"),
                      scale = ~ 1, dependence = ~ 1, root = 1, shape = ~ 1,
                      weights = NULL, missing = c("fail", "likelihood", "complete"),
                      split = NULL, seed = NULL, copula = "pcc",
                      covariates = c("none", "backward"), ...) {
  criterion <- match.arg(criterion)
  missing <- match.arg(missing)
  if (any(c("truncation", "censoring") %in% names(list(...)))) {
    stop("tc_select() does not support truncation or censoring.")
  }
  covariates <- match.arg(covariates)
  if (identical(copula, "auto")) copula <- c("gaussian", "student", "pcc")
  if (!is.character(copula) || !length(copula) || anyDuplicated(copula) ||
      any(!copula %in% c("gaussian", "student", "pcc"))) {
    stop("copula must be \"auto\" or among \"gaussian\", \"student\", \"pcc\".")
  }
  if (is.list(margins)) {
    if (length(margins) != 3L || !all(vapply(margins, is.character, logical(1))) ||
        any(lengths(margins) == 0L) || any(!unlist(margins) %in% .tc_margin_names)) {
      stop("margins must be three margin names or a list of three vectors of candidate names.")
    }
  }
  margins0 <- margins; formulas0 <- formulas
  sel_margins <- is.list(margins) && any(lengths(margins) > 1L)
  all_cand <- identical(candidates, "all")
  if (all_cand) candidates <- NULL
  if (!is.null(candidates) && any(!candidates %in% .tc_pair_families)) {
    stop("candidates must be \"all\" or among: ", paste(.tc_pair_families, collapse = ", "))
  }
  if (length(root) != 1L || !(root %in% 1:3 || identical(root, "auto"))) {
    stop("root must be 1, 2, 3 or \"auto\".")
  }
  # settings kept with the fit so that tc_bootstrap(select = TRUE) can repeat
  # the selection (0.14.0)
  settings <- list(candidates = if (all_cand) "all" else candidates, criterion = criterion,
                   root = root, dots = list(...), split = split,
                   margins = margins0, formulas = formulas0, copula = copula,
                   covariates = covariates, scale = scale, dependence = dependence,
                   shape = shape)
  if (!is.null(split)) {
    # sample splitting (0.14.0): select on one part, estimate on the other
    if (!is.numeric(split) || length(split) != 1L || !is.finite(split) ||
        split <= 0 || split >= 1) {
      stop("split must be a single number strictly between 0 and 1.")
    }
    if (!is.null(weights) && length(weights) != nrow(data)) {
      stop("weights must be a vector of non-negative numbers, one per row of data.")
    }
    if (!is.null(seed)) set.seed(seed)
    n0 <- nrow(data)
    in_sel <- seq_len(n0) %in% sample.int(n0, round(split * n0))
    if (sum(in_sel) < 20L || sum(!in_sel) < 20L) {
      stop("split leaves fewer than 20 rows in one of the two parts.")
    }
    sel_args <- list(...)
    sel_args$hessian <- FALSE
    sel <- do.call(tc_select, c(list(
      formulas, data[in_sel, , drop = FALSE], margins = margins,
      candidates = settings$candidates, criterion = criterion, scale = scale,
      dependence = dependence, root = root, shape = shape, weights = weights[in_sel],
      missing = missing, copula = copula, covariates = covariates), sel_args))
    fit <- tricopula(sel$formulas, data[!in_sel, , drop = FALSE], margins = sel$margins,
                     copula = sel$copula, families = sel$families, scale = scale,
                     dependence = dependence,
                     root = if (is.null(sel$root)) 1 else sel$root, shape = shape,
                     weights = weights[!in_sel], missing = missing, ...)
    fit$selection <- settings
    out <- sel
    out$fit <- fit
    out$split <- split
    out$rows_select <- which(in_sel); out$rows_estimate <- which(!in_sel)
    out$selection_fit <- sel$fit
    return(out)
  }
  # margins and covariates (0.14.1), by the marginal fits
  pre <- NULL
  if (sel_margins || covariates == "backward") {
    fl <- function(f) if (inherits(f, "formula")) list(f) else f
    pre <- .tc_select_margins(formulas, data,
                              if (is.list(margins)) margins else as.list(margins),
                              covariates, criterion, scale, shape, weights, missing,
                              c(fl(dependence), fl(list(...)$df)))
    margins <- pre$margins; formulas <- pre$formulas
  } else if (is.list(margins)) margins <- unlist(margins)
  finish <- function(out) {
    out$copula <- "pcc"
    ell <- setdiff(copula, "pcc")
    if (length(ell)) {
      # copula type (0.14.1): joint fits compared by the criterion
      crit_of <- function(f) {
        if (is.null(f) || inherits(f, "try-error")) return(Inf)
        if (criterion == "AIC") stats::AIC(f) else stats::BIC(f)
      }
      fits <- lapply(ell, function(cc) {
        try(tricopula(formulas, data, margins = margins, copula = cc, scale = scale,
                      dependence = dependence, shape = shape, weights = weights,
                      missing = missing, ...), silent = TRUE)
      })
      names(fits) <- ell
      if (!is.null(out$fit)) fits <- c(fits, list(pcc = out$fit))
      cr <- vapply(fits, crit_of, numeric(1))
      if (!any(is.finite(cr))) stop("None of the copula types could be fitted.")
      best <- names(fits)[which.min(cr)]
      ctab <- data.frame(
        copula = names(fits),
        logLik = vapply(fits, function(f) if (inherits(f, "try-error")) NA_real_ else f$logLik,
                        numeric(1)),
        npar = vapply(fits, function(f) if (inherits(f, "try-error")) NA_integer_ else
          f$parameter_count, integer(1)),
        criterion = cr, selected = names(fits) == best, row.names = NULL)
      names(ctab)[names(ctab) == "criterion"] <- criterion
      out$copulas <- ctab
      out$copula <- best
      if (best != "pcc") {
        out$pcc_families <- out$families
        out$fit <- fits[[best]]
        out["families"] <- list(NULL); out["root"] <- list(NULL)
      }
    }
    out$fit$selection <- settings
    out$margins <- margins; out$formulas <- formulas
    out$margin_table <- pre$margin_table; out$covariate_table <- pre$covariate_table
    out
  }
  if (!"pcc" %in% copula) {
    return(finish(list(families = NULL, table = NULL, fit = NULL, criterion = criterion,
                       root = NULL)))
  }
  # marginal fits and their PITs (the first stage of the two-stage estimator)
  base <- tricopula(formulas, data, margins = margins, copula = "gaussian",
                    hessian = FALSE, polish = FALSE, scale = scale, shape = shape,
                    control = list(maxit = 0),   # only the marginal stage is used

                    weights = weights, missing = missing)
  idx <- base$index
  U0 <- sapply(1:3, function(j) {
    X <- base$design[[j]]; Z <- base$scale_design[[j]]
    b <- base$ifm$coefficients[idx$margin[[j]]]
    if (length(idx$shapel[[j]])) {
      pe <- ncol(X) + ncol(Z)
      return(.tc_margin_eval(base$response[[j]], X, margins[j], b[seq_len(ncol(X))],
                             drop(Z %*% b[ncol(X) + seq_len(ncol(Z))]),
                             shape = drop(base$shape_design[[j]] %*% b[-seq_len(pe)]))$u)
    }
    .tc_margin_eval(base$response[[j]], X, margins[j], b[seq_len(ncol(X))],
                    drop(Z %*% b[-seq_len(ncol(X))]))$u
  })
  if (!is.null(base$observed)) U0[!base$observed] <- NA_real_
  n <- nrow(U0)
  pen <- if (criterion == "AIC") 2 else log(n)
  pick <- function(u, v, label) {
    # missing responses: the pairs observed together (unweighted)
    cc <- !is.na(u) & !is.na(v)
    u <- u[cc]; v <- v[cc]
    cand <- candidates
    if (is.null(cand)) {
      cand <- if (stats::cor(u, v, method = "kendall") >= 0) {
        c("gaussian", "student", "clayton", "gumbel", "clayton180", "gumbel180",
          if (all_cand) c("frank", "joe", "joe180"))
      } else c("gaussian", "student", "clayton90", "gumbel90", "clayton270", "gumbel270",
               if (all_cand) c("frank", "joe90", "joe270"))
    }
    fits <- lapply(cand, function(f) .tc_fit_pair(u, v, f))
    tab <- data.frame(pair = label, family = cand,
                      logLik = vapply(fits, `[[`, numeric(1), "logLik"),
                      npar = vapply(fits, `[[`, numeric(1), "npar"))
    tab$criterion <- -2 * tab$logLik + pen * tab$npar
    best <- which.min(tab$criterion)
    tab$selected <- seq_len(nrow(tab)) == best
    list(tab = tab, family = cand[best], fit = fits[[best]])
  }
  one_root <- function(r) {
    U <- U0[, .tc_pcc_order(r), drop = FALSE]
    pl <- .tc_pcc_labels(r)$short
    s12 <- pick(U[, 1], U[, 2], pl[1])
    s13 <- pick(U[, 1], U[, 3], pl[2])
    a <- .tc_pair(U[, 1], U[, 2], s12$fit$eta, s12$family, s12$fit$tau)$h
    b <- .tc_pair(U[, 1], U[, 3], s13$fit$eta, s13$family, s13$fit$tau)$h
    s23 <- pick(a, b, pl[3])
    fam <- c(s12$family, s13$family, s23$family)
    tab <- rbind(s12$tab, s13$tab, s23$tab)
    names(tab)[names(tab) == "criterion"] <- criterion
    fit <- tricopula(formulas, data, margins = margins, copula = "pcc",
                     families = fam, scale = scale, dependence = dependence,
                     root = r, shape = shape, weights = weights, missing = missing, ...)
    fit$selection <- settings
    list(families = fam, table = tab, fit = fit, criterion = criterion,
         root = as.integer(r))
  }
  if (!identical(root, "auto")) return(finish(one_root(root)))
  res <- lapply(1:3, one_root)
  crit <- vapply(res, function(z) {
    if (criterion == "AIC") stats::AIC(z$fit) else stats::BIC(z$fit)
  }, numeric(1))
  best <- which.min(crit)
  out <- res[[best]]
  out$roots <- data.frame(
    root = 1:3,
    families = vapply(res, function(z) paste(z$families, collapse = ","), character(1)),
    logLik = vapply(res, function(z) z$fit$logLik, numeric(1)),
    npar = vapply(res, function(z) z$fit$parameter_count, integer(1)),
    criterion = crit, selected = 1:3 == best)
  names(out$roots)[names(out$roots) == "criterion"] <- criterion
  finish(out)
}

# Selection of the margin families and of the terms of the mean formulas by
# marginal maximum likelihood (0.14.1). `cand` is a list of three vectors of
# candidate margins; `extra` holds further formulas whose variables decide
# which rows tricopula() would drop.
.tc_select_margins <- function(formulas, data, cand, covariates, criterion, scale, shape,
                               weights, missing, extra = NULL) {
  fl <- function(f) if (inherits(f, "formula")) list(f) else f
  if (any(vapply(c(formulas, fl(scale), fl(shape)), .tc_has_smooth, logical(1)))) {
    stop("Selection of margins or covariates is not available for formulas with smooth terms.")
  }
  wt <- if (is.null(weights)) NULL else as.numeric(weights)
  if (!is.null(wt) && length(wt) != nrow(data)) {
    stop("weights must be a vector of non-negative numbers, one per row of data.")
  }
  if (missing != "fail") {
    # the rows tricopula() uses
    vars_y <- lapply(formulas, function(f) intersect(all.vars(f[[2]]), names(data)))
    vars_all <- unique(unlist(lapply(c(formulas, fl(scale), fl(shape), extra), all.vars)))
    vars_x <- setdiff(intersect(vars_all, names(data)), unlist(vars_y))
    keep <- .tc_missing_rows(data, vars_y, vars_x, missing)
    data <- data[keep, , drop = FALSE]
    if (!is.null(wt)) wt <- wt[keep]
  }
  n <- nrow(data)
  zz <- .tc_make_W(scale, data, n, what = "scale")
  qq <- .tc_make_W(shape, data, n, what = "shape")
  mtab <- list(); ctab <- list()
  out_m <- character(3); out_f <- formulas
  for (j in 1:3) {
    Z <- zz[[j]]$W; Q <- qq[[j]]$W
    fitj <- function(f, m) {
      mm <- .tc_make_X(f, data, na_ok = missing == "likelihood")
      y <- mm$y; ob <- !is.na(y)
      if (!ncol(mm$X)) return(NULL)
      ok <- try(.tc_check_support(y[ob], m), silent = TRUE)
      if (inherits(ok, "try-error")) return(NULL)
      if (!all(ob)) y[!ob] <- stats::median(y[ob])
      w <- if (is.null(wt) && all(ob)) NULL else (if (is.null(wt)) 1 else wt) * ob
      r <- try(.tc_fit_margin(y, mm$X, m, Z, if (m %in% .tc_margin_names3) Q, w),
               silent = TRUE)
      if (inherits(r, "try-error") || !is.finite(r$logLik)) return(NULL)
      pen <- if (criterion == "AIC") 2 else log(sum(ob))
      list(logLik = r$logLik, npar = length(r$par),
           crit = -2 * r$logLik + pen * length(r$par))
    }
    cj <- cand[[j]]
    fits <- lapply(cj, function(m) fitj(formulas[[j]], m))
    crit <- vapply(fits, function(z) if (is.null(z)) Inf else z$crit, numeric(1))
    if (!any(is.finite(crit))) stop("No candidate margin could be fitted for response ", j, ".")
    best <- which.min(crit)
    out_m[j] <- cj[best]
    if (length(cj) > 1L) {
      mtab[[j]] <- data.frame(
        margin = j, family = cj,
        logLik = vapply(fits, function(z) if (is.null(z)) NA_real_ else z$logLik, numeric(1)),
        npar = vapply(fits, function(z) if (is.null(z)) NA_integer_ else z$npar, integer(1)),
        criterion = crit, selected = seq_along(cj) == best)
    }
    if (covariates == "backward") {
      f <- formulas[[j]]
      tt <- stats::terms(f)
      all_terms <- attr(tt, "term.labels")
      ft <- function(tl) {
        g <- stats::reformulate(if (length(tl)) tl else "1", response = f[[2]],
                                intercept = attr(tt, "intercept") == 1)
        environment(g) <- environment(f)
        g
      }
      keep_t <- all_terms
      cur <- crit[best]
      while (length(keep_t)) {
        scope <- stats::drop.scope(stats::terms(ft(keep_t)))
        if (!length(scope)) break
        cr <- vapply(scope, function(s) {
          z <- fitj(ft(setdiff(keep_t, s)), out_m[j])
          if (is.null(z)) Inf else z$crit
        }, numeric(1))
        if (min(cr) >= cur) break
        keep_t <- setdiff(keep_t, scope[which.min(cr)])
        cur <- min(cr)
      }
      if (length(keep_t) < length(all_terms)) out_f[[j]] <- ft(keep_t)
      if (length(all_terms)) {
        ctab[[j]] <- data.frame(margin = j, term = all_terms,
                                selected = all_terms %in% keep_t)
      }
    }
  }
  mt <- if (length(mtab)) do.call(rbind, mtab)
  if (!is.null(mt)) {
    names(mt)[names(mt) == "criterion"] <- criterion
    rownames(mt) <- NULL
  }
  list(margins = out_m, formulas = out_f, margin_table = mt,
       covariate_table = if (length(ctab)) do.call(rbind, ctab))
}
