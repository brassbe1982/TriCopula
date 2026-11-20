# Non-simplified vine (0.15.0) ----------------------------------------------
#
# In a simplified C-vine the conditional pair copula c_ab|r does not depend on
# the value of the conditioning (root) variable. With simplified = FALSE its
# predictor becomes
#     eta_ab|r,i = w_i' delta + f(u_r,i),     u_r = F_r(y_r | x),
# with f a penalized cubic B-spline on (0, 1) (a P-spline: k equally spaced
# basis functions, second-order difference penalty) under the constraint
# int_0^1 f(u) du = 0, so that the intercept keeps its meaning as the average
# predictor over the conditioning value. The knots are fixed, not taken from
# the data. In joint estimation u_r is itself a function of the parameters of
# margin r; the score therefore has the extra term
#     d logc / d u_r += (d logc / d eta_ab|r) * f'(u_r).
# The coefficients of f are the last elements of the parameter vector.

.tc_ns_setup <- function(k = 10L) {
  k <- suppressWarnings(as.integer(k))
  if (length(k) != 1L || is.na(k) || k < 4L || k > 40L) {
    stop("k_cond must be a whole number between 4 and 40.", call. = FALSE)
  }
  nseg <- k - 3L
  knots <- seq(-3L, nseg + 3L) / nseg
  full <- function(u, d = 0L) {
    splines::splineDesign(knots, u, ord = 4L, derivs = rep(d, length(u)), outer.ok = TRUE)
  }
  # int_0^1 B_j(u) du: Gauss-Legendre on each knot interval (exact for cubics)
  gl <- .tc_gl(16L)
  nodes <- as.vector(outer(gl$x, seq_len(nseg) - 1L, "+")) / nseg
  cvec <- colSums(full(nodes) * (rep(gl$w, nseg) / nseg))
  Zc <- qr.Q(qr(matrix(cvec, ncol = 1L)), complete = TRUE)[, -1L, drop = FALSE]
  D <- diff(diag(k), differences = 2L)
  S <- crossprod(D %*% Zc)
  S <- (S + t(S)) / 2
  S <- S / max(abs(S))
  list(k = k, ncoef = k - 1L, knots = knots, Zc = Zc, S = S,
       basis = function(u, deriv = TRUE) {
         u <- pmin(pmax(u, 0), 1)
         list(B = full(u) %*% Zc, dB = if (deriv) full(u, 1L) %*% Zc)
       })
}

#' Test of the simplifying assumption of a pair-copula construction
#'
#' Likelihood-ratio test of a simplified C-vine against the vine in which the
#' predictor of the conditional pair copula is a function of the
#' conditioning value (see `simplified` in [tricopula()]).
#'
#' @details
#' The alternative is fitted without penalty (smoothing parameter 0) with a
#' small spline basis, so that the two models are nested parametric models
#' and the statistic has a chi-square reference distribution with `k - 1`
#' degrees of freedom. The test has power against alternatives that a spline
#' with `k - 1` coefficients can approximate; `k = 4` detects monotone and
#' simple curved departures.
#'
#' @param object A simplified `tricopula` fit with `copula = "pcc"`.
#' @param k Number of B-spline basis functions of the alternative (at least 4).
#' @return A data frame with the two log-likelihoods, the statistic, its
#'   degrees of freedom and the p-value; attribute `"fit"` holds the fit of
#'   the alternative (without covariance matrix).
#' @export
tc_simplifying_test <- function(object, k = 4L) {
  .tc_check(object)
  if (object$copula != "pcc") stop("The test needs a fit with copula = \"pcc\".")
  if (isFALSE(object$simplified)) stop("The test needs the simplified fit.")
  ns <- .tc_ns_setup(k)
  full <- .tc_refit(object, object$data, sp = c(unname(object$sp), 0), hessian = FALSE,
                    start = c(object$coefficients, rep(0, ns$ncoef)),
                    simplified = FALSE, k_cond = k)
  if (inherits(full, "try-error")) stop("The non-simplified fit failed: ", full)
  stat <- max(2 * (full$logLik - object$logLik), 0)
  out <- data.frame(logLik_simplified = object$logLik, logLik_nonsimplified = full$logLik,
                    statistic = stat, df = ns$ncoef,
                    p.value = stats::pchisq(stat, ns$ncoef, lower.tail = FALSE))
  attr(out, "fit") <- full
  out
}
