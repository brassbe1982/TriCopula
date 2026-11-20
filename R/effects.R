# Marginal effects ---------------------------------------------------------

# Quantity of interest per covariate row for a given coefficient vector:
# an n x 3 matrix of conditional means, or an n x 1 matrix of joint
# probabilities P(Y1 <= t1, Y2 <= t2, Y3 <= t3 | x).
.tc_quantity <- function(par, object, newdata, type, thresholds) {
  object$coefficients[] <- par
  if (object$copula == "student" && !isFALSE(object$df_constant)) {
    object$df <- 2.05 + exp(par[[object$index$df]])
  }
  if (type == "mean") {
    as.matrix(predict.tricopula(object, newdata, type = "mean"))
  } else {
    matrix(tc_jointprob(object, thresholds, newdata), ncol = 1,
           dimnames = list(NULL, "jointprob"))
  }
}

#' Average marginal effects
#'
#' Average effect of one covariate on the three conditional means
#' `E(Y_j | x)` or on the joint probability
#' `P(Y1 <= t1, Y2 <= t2, Y3 <= t3 | x)`, with delta-method standard errors.
#' The covariate acts through every part of the model in which it appears
#' (means, scale/shape and dependence predictors).
#'
#' @details
#' For a numeric covariate the effect is the average derivative, computed by
#' a central difference of width `2 * h` around each observed value. If
#' `values = c(a, b)` is given, the effect is instead the average change
#' when the covariate moves from `a` to `b` for every row (use this for
#' factors, dummies and discrete changes). Standard errors treat the
#' covariate rows as fixed and propagate the uncertainty of the estimated
#' coefficients through a numerical Jacobian.
#'
#' `type = "jointprob"` evaluates a trivariate integral per distinct row and
#' per coefficient, so it is slow; pass a small `newdata`.
#'
#' @param object A `tricopula` fit.
#' @param variable Name of the covariate.
#' @param newdata Data frame of covariate rows to average over; defaults to
#'   the estimation data.
#' @param type `"mean"` or `"jointprob"`.
#' @param thresholds Three thresholds, required for `type = "jointprob"`.
#' @param values Optional vector of two values for a discrete change.
#' @param h Half-width of the finite difference for a numeric covariate;
#'   defaults to `1e-4` times its standard deviation.
#' @param level Confidence level.
#' @param vcov_type `"hessian"` or `"sandwich"`; default as for `vcov()`.
#' @return A data frame with one row per outcome: the effect, its standard
#'   error, a Wald interval and p-value.
#' @export
tc_effects <- function(object, variable, newdata = NULL,
                       type = c("mean", "jointprob"), thresholds = NULL,
                       values = NULL, h = NULL, level = 0.95,
                       vcov_type = NULL) {
  .tc_check(object)
  type <- match.arg(type)
  vcov_type <- .tc_vcov_type(object, vcov_type)
  if (is.null(newdata)) newdata <- object$data
  if (is.null(newdata)) stop("Supply newdata: the fit does not contain the data.")
  if (!is.character(variable) || length(variable) != 1L ||
      !variable %in% names(newdata)) {
    stop("variable must be the name of one column of the data.")
  }
  if (type == "jointprob" && (!is.numeric(thresholds) || length(thresholds) != 3L)) {
    stop("thresholds must be a numeric vector of length 3 for type = \"jointprob\".")
  }
  x <- newdata[[variable]]
  lo <- hi <- newdata
  if (!is.null(values)) {
    if (length(values) != 2L) stop("values must have length 2.")
    cast <- function(v) if (is.factor(x)) factor(v, levels = levels(x)) else v
    lo[[variable]] <- cast(rep(values[1], nrow(newdata)))
    hi[[variable]] <- cast(rep(values[2], nrow(newdata)))
    if (anyNA(lo[[variable]]) || anyNA(hi[[variable]])) stop("values are not valid levels of the variable.")
    denom <- 1
  } else {
    if (!is.numeric(x)) stop("For a non-numeric variable supply values = c(from, to).")
    if (is.null(h)) h <- 1e-4 * max(stats::sd(x), 1e-8, na.rm = TRUE)
    if (!is.finite(h) || h <= 0) h <- 1e-4
    lo[[variable]] <- x - h
    hi[[variable]] <- x + h
    denom <- 2 * h
  }
  ame <- function(par) {
    colMeans(.tc_quantity(par, object, hi, type, thresholds) -
               .tc_quantity(par, object, lo, type, thresholds)) / denom
  }
  est <- ame(object$coefficients)
  se <- rep(NA_real_, length(est))
  if (!is.null(object$vcov)) {
    V <- stats::vcov(object, type = vcov_type)
    V[is.na(V)] <- 0
    J <- if (type == "mean") {
      numDeriv::jacobian(ame, object$coefficients)
    } else {
      numDeriv::jacobian(ame, object$coefficients, method = "simple",
                         method.args = list(eps = 1e-4))
    }
    se <- sqrt(pmax(diag(J %*% V %*% t(J)), 0))
  }
  zq <- stats::qnorm(1 - (1 - level) / 2)
  out <- data.frame(
    outcome = if (type == "mean") object$response_names else
      sprintf("P(%s <= %s)", paste(object$response_names, collapse = ", "),
              paste(format(thresholds), collapse = ", ")),
    effect = unname(est), std.error = se,
    lower = unname(est) - zq * se, upper = unname(est) + zq * se,
    p.value = 2 * stats::pnorm(abs(unname(est) / se), lower.tail = FALSE))
  attr(out, "variable") <- variable
  attr(out, "kind") <- if (is.null(values)) "average derivative" else
    paste0("average change from ", values[1], " to ", values[2])
  out
}
