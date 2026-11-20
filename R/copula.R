# Trivariate elliptical copulas ------------------------------------------

# Parameter layout:
#   margin 1: beta_1, gamma_1 | margin 2: beta_2, gamma_2 | margin 3: ...
#   (gamma_j: coefficients of the log-scale / log-shape predictor, e_j of them)
#   delta_12 (q1), delta_13 (q2), delta_23|1 (q3) [, log(df - 2.05)]
# where atanh(rho12_i) = w1_i' delta_12 etc. With intercept-only dependence
# (q = 1, 1, 1) this is exactly the 0.1.0/0.2.0 layout.
#   r: number of coefficients of each log(df - 2.05) predictor (0.10.0); the
#   Student-t copula uses r[1], a PCC r[k] for each Student-t pair k.
#   s: number of coefficients of the third (shape) predictor of each margin
#   (0.11.0); 0 for two-parameter margins. They follow gamma_j inside the
#   block of margin j, so the layout of two-parameter models is unchanged.
.tc_index <- function(p, copula, q = c(1L, 1L, 1L), e = c(1L, 1L, 1L),
                      r = c(1L, 1L, 1L), s = c(0L, 0L, 0L), ns = 0L) {
  idx <- list(beta = vector("list", 3L), extral = vector("list", 3L),
              shapel = vector("list", 3L))
  pos <- 0L
  for (j in 1:3) {
    idx$beta[[j]] <- pos + seq_len(p[j])
    pos <- pos + p[j]
    idx$extral[[j]] <- pos + seq_len(e[j])
    pos <- pos + e[j]
    idx$shapel[[j]] <- pos + seq_len(s[j])
    pos <- pos + s[j]
  }
  idx$e <- as.integer(e)
  idx$depl <- vector("list", 3L)
  for (k in 1:3) {
    idx$depl[[k]] <- pos + seq_len(q[k])
    pos <- pos + q[k]
  }
  idx$dep <- unlist(idx$depl)
  idx$q <- as.integer(q)
  # one df parameter for the Student-t copula; one per Student-t pair in a PCC
  nd <- if (copula == "student") as.integer(r[1]) else if (copula == "pcc") {
    ifelse(.tc_pair_base(attr(copula, "families")) == "student", as.integer(r), 0L)
  } else integer(0)
  dend <- pos + cumsum(nd)
  idx$dfl <- lapply(seq_along(nd), function(k) dend[k] - nd[k] + seq_len(nd[k]))
  idx$df <- pos + seq_len(sum(nd))
  idx$npar <- pos + length(idx$df)
  # coefficients of the spline in the conditioning value of a non-simplified
  # vine (0.15.0): the last elements, so that every other position is unchanged
  idx$ns <- idx$npar + seq_len(ns)
  idx$npar <- idx$npar + as.integer(ns)
  idx$margin <- lapply(1:3, function(j) c(idx$beta[[j]], idx$extral[[j]], idx$shapel[[j]]))
  idx
}

# Map unconstrained dependence predictors (scalars or vectors) to
# correlations through the partial correlation rho23|1 (C-vine order 1,2,3):
#   rho23 = rho12 rho13 + sqrt((1 - rho12^2)(1 - rho13^2)) rho23|1.
# det R = (1 - rho12^2)(1 - rho13^2)(1 - rho23|1^2) > 0 for all real inputs.
.tc_dep_vec <- function(a, b, c) {
  r12 <- tanh(a)
  r13 <- tanh(b)
  p23 <- tanh(c)
  s <- sqrt((1 - r12^2) * (1 - r13^2))
  list(r12 = r12, r13 = r13, r23 = r12 * r13 + s * p23, p23 = p23, s = s)
}

# Scalar version returning the correlation matrix (constant dependence).
.tc_dep <- function(theta, copula) {
  d <- .tc_dep_vec(theta[[1]], theta[[2]], theta[[3]])
  d$R <- .tc_R(d$r12, d$r13, d$r23)
  d$df <- if (copula == "student") 2.05 + exp(theta[[4]]) else Inf
  d
}

.tc_R <- function(r12, r13, r23) {
  matrix(c(1, r12, r13, r12, 1, r23, r13, r23, 1), 3, 3)
}

# Largest admissible |atanh(rho)| (see .tc_cop_eval).
.tc_eta_max <- 6

# Inverse map: correlation matrix -> unconstrained parameters.
.tc_dep_inverse <- function(R, df = NULL) {
  lim <- function(r) pmin(pmax(r, -0.999), 0.999)
  r12 <- lim(R[1, 2]); r13 <- lim(R[1, 3])
  p23 <- lim((R[2, 3] - r12 * r13) / sqrt((1 - r12^2) * (1 - r13^2)))
  th <- c(atanh(r12), atanh(r13), atanh(p23))
  if (!is.null(df) && is.finite(df)) th <- c(th, log(df - 2.05))
  th
}

# Log copula density for an n x 3 matrix of uniforms.
#   Gaussian: -1/2 log|R| - 1/2 z'(R^-1 - I) z,            z = qnorm(u)
#   Student : log t_3(z; R, df) - sum_j log t_1(z_j; df),  z = qt(u, df)
# R is a 3 x 3 correlation matrix or a list with elements r12, r13, r23,
# each a scalar or a vector of length n (observation-specific dependence).
# The 3 x 3 inverse and determinant are written out in closed form so that
# everything is vectorized over observations.
# With deriv = TRUE also returns d logc / d u (n x 3) and
# d logc / d (rho12, rho13, rho23) (n x 3).
.tc_logcop <- function(U, R, df, copula, deriv = FALSE) {
  if (is.matrix(R)) R <- list(r12 = R[1, 2], r13 = R[1, 3], r23 = R[2, 3])
  r12 <- R$r12; r13 <- R$r13; r23 <- R$r23
  # det R = (1 - rho12^2)(1 - rho13^2)(1 - rho23|1^2): the product form avoids
  # the cancellation of 1 + 2 r12 r13 r23 - r12^2 - r13^2 - r23^2 near
  # singular matrices.
  p23 <- if (!is.null(R$p23)) R$p23 else
    (r23 - r12 * r13) / sqrt((1 - r12^2) * (1 - r13^2))
  det <- (1 - r12^2) * (1 - r13^2) * (1 - p23^2)
  i11 <- (1 - r23^2) / det; i22 <- (1 - r13^2) / det; i33 <- (1 - r12^2) / det
  i12 <- (r13 * r23 - r12) / det
  i13 <- (r12 * r23 - r13) / det
  i23 <- (r12 * r13 - r23) / det
  ld <- log(det)
  Z <- if (copula == "gaussian") stats::qnorm(U) else stats::qt(U, df = df)
  B <- cbind(Z[, 1] * i11 + Z[, 2] * i12 + Z[, 3] * i13,
             Z[, 1] * i12 + Z[, 2] * i22 + Z[, 3] * i23,
             Z[, 1] * i13 + Z[, 2] * i23 + Z[, 3] * i33)
  if (copula == "gaussian") {
    lc <- -0.5 * ld - 0.5 * rowSums((B - Z) * Z)
    w <- 1
    if (deriv) du <- -(B - Z) / stats::dnorm(Z)
  } else {
    q <- rowSums(B * Z)
    log_uni <- stats::dt(Z, df = df, log = TRUE)
    lc <- lgamma((df + 3) / 2) - lgamma(df / 2) - 1.5 * log(df * pi) -
      0.5 * ld - ((df + 3) / 2) * log1p(q / df) - rowSums(log_uni)
    w <- (df + 3) / (df + q)
    if (deriv) du <- (-w * B + (df + 1) * Z / (df + Z^2)) / exp(log_uni)
  }
  out <- list(lc = lc)
  if (deriv) {
    out$du <- du
    out$dR <- cbind(-i12 + w * B[, 1] * B[, 2],
                    -i13 + w * B[, 1] * B[, 3],
                    -i23 + w * B[, 2] * B[, 3])
  }
  out
}

# Copula log-density as a function of the dependence coefficient vector
# theta = (delta_12, delta_13, delta_23|1 [, log(df - 2.05)]) for given
# uniforms U and dependence design matrices W (list of three, n rows each).
# With deriv = TRUE also the n x length(theta) score matrix and d logc / d u.
# V (0.10.0): optional list of design matrices of the log(df - 2.05)
# predictors; NULL or an intercept column gives the constant-df code path of
# earlier versions. For a PCC, attribute `root` of `copula` selects the vine
# root: the columns of U are passed to the vine in the order (root, others).
.tc_cop_eval <- function(theta, U, W, copula, deriv = FALSE, V = NULL) {
  # missing responses (0.13.0): observed-data copula terms, see missing.R
  if (anyNA(U)) return(.tc_cop_eval_missing(theta, U, W, copula, deriv = deriv, V = V))
  if (copula == "pcc") {
    root <- attr(copula, "root")
    if (is.null(root) || root == 1L) {
      return(.tc_pcc_eval(theta, U, W, attr(copula, "families"), deriv = deriv, V = V,
                          ns = attr(copula, "ns")))
    }
    ord <- .tc_pcc_order(root)
    out <- .tc_pcc_eval(theta, U[, ord, drop = FALSE], W, attr(copula, "families"),
                        deriv = deriv, V = V, ns = attr(copula, "ns"))
    if (deriv) {
      du <- out$du
      out$du[, ord] <- du
    }
    return(out)
  }
  q <- vapply(W, ncol, integer(1))
  ends <- cumsum(q)
  pos <- lapply(1:3, function(k) ends[k] - q[k] + seq_len(q[k]))
  eta <- lapply(1:3, function(k) drop(W[[k]] %*% theta[pos[[k]]]))
  dep <- .tc_dep_vec(eta[[1]], eta[[2]], eta[[3]])
  dfc <- copula != "student" || .tc_df_const(V[[1]])
  dpos <- if (copula == "student") ends[3] + seq_len(if (dfc) 1L else ncol(V[[1]]))
  tau <- if (copula != "student") NULL else if (dfc) theta[[dpos]] else
    as.vector(V[[1]] %*% theta[dpos])
  df <- if (copula == "student") 2.05 + exp(tau) else Inf
  # Numerical domain: |rho| <= tanh(6) = 0.999988. Beyond it the correlation
  # matrix is singular in floating point and the density evaluates to
  # rounding noise (which can look like a very large likelihood), so such
  # points are rejected.
  if (any(abs(unlist(eta)) > .tc_eta_max)) {
    return(list(lc = rep(-Inf, nrow(U)), dep = dep, df = df))
  }
  cp <- .tc_logcop(U, dep, df, copula, deriv = deriv)
  out <- list(lc = cp$lc, dep = dep, df = df)
  if (deriv) {
    r12 <- dep$r12; r13 <- dep$r13; p23 <- dep$p23
    ds12 <- -r12 * sqrt((1 - r13^2) / (1 - r12^2))
    ds13 <- -r13 * sqrt((1 - r12^2) / (1 - r13^2))
    S <- matrix(0, nrow(U), length(theta))
    S[, pos[[1]]] <- W[[1]] *
      ((cp$dR[, 1] + cp$dR[, 3] * (r13 + p23 * ds12)) * (1 - r12^2))
    S[, pos[[2]]] <- W[[2]] *
      ((cp$dR[, 2] + cp$dR[, 3] * (r12 + p23 * ds13)) * (1 - r13^2))
    S[, pos[[3]]] <- W[[3]] * (cp$dR[, 3] * dep$s * (1 - p23^2))
    if (copula == "student") {
      # d z / d df has no closed form: central difference in log(df - 2.05).
      h <- 1e-5
      lp <- .tc_logcop(U, dep, 2.05 + exp(tau + h), copula)$lc
      lm <- .tc_logcop(U, dep, 2.05 + exp(tau - h), copula)$lc
      S[, dpos] <- if (dfc) (lp - lm) / (2 * h) else V[[1]] * ((lp - lm) / (2 * h))
    }
    out$score <- S
    out$du <- cp$du
  }
  out
}

# Joint log-likelihood; with deriv = TRUE also the n x npar matrix of
# per-observation scores with respect to the unconstrained parameters.
# Q (0.11.0): list of design matrices of the third (shape) predictor of
# three-parameter margins; unused entries are ignored.
# wt, obs (0.13.0): prior weights (one per row) and the n x 3 logical matrix
# of observed responses. A missing response (its y holds a placeholder)
# contributes no density, has PIT NA and zero derivatives. With wt the
# log-likelihood is sum(wt * ll_i) and each score row is multiplied by its
# weight; ll_i itself stays unweighted.
.tc_eval <- function(par, y, X, margins, copula, idx, deriv = FALSE, W, Z,
                     V = NULL, Q = NULL, wt = NULL, obs = NULL) {
  # truncation and censoring (0.16.0): carried as attribute "ct" of y, so
  # that every caller of .tc_eval() uses them; see censoring.R
  ct <- attr(y, "ct")
  mc <- NULL
  if (is.null(ct)) {
    me <- lapply(1:3, function(j) {
      is3 <- idx$shapel[[j]]
      .tc_margin_eval(y[[j]], X[[j]], margins[j], par[idx$beta[[j]]],
                      drop(Z[[j]] %*% par[idx$extral[[j]]]), deriv = deriv,
                      shape = if (length(is3)) drop(Q[[j]] %*% par[is3]))
    })
  } else {
    mc <- .tc_margins_ct(par, y, X, margins, idx, Z, Q, deriv, ct)
    me <- mc$me
    if (is.null(mc$cen)) mc <- NULL
  }
  if (!is.null(obs)) {
    for (j in 1:3) {
      mis <- !obs[, j]
      if (!any(mis)) next
      me[[j]]$logf[mis] <- 0
      me[[j]]$u[mis] <- NA_real_
      for (nm in intersect(c("dl_eta", "dl_extra", "dl_shape", "du_eta", "du_extra",
                             "du_shape"), names(me[[j]]))) {
        me[[j]][[nm]][mis] <- 0
      }
    }
  }
  wsum <- function(x) if (is.null(wt)) sum(x) else sum(wt * x)
  U <- do.call(cbind, lapply(me, `[[`, "u"))
  cp <- if (is.null(mc)) {
    .tc_cop_eval(par[c(idx$dep, idx$df, idx$ns)], U, W, copula, deriv = deriv, V = V)
  } else {
    .tc_cop_eval_cens(par[c(idx$dep, idx$df, idx$ns)], mc$Lo, mc$Up, mc$cen, W, copula,
                      deriv = deriv, V = V, nodes = ct$nodes,
                      method = if (is.null(ct$method)) "quadrature" else ct$method)
  }
  ll_i <- rowSums(do.call(cbind, lapply(me, `[[`, "logf"))) + cp$lc
  if (any(!is.finite(ll_i))) return(list(loglik = -Inf, ll_i = ll_i))
  dep <- cp$dep
  if (copula == "pcc") {
    fam <- dep$families
    nat <- sapply(1:3, function(k) .tc_pair_natural(dep$eta[[k]], fam[k]))
    kt <- sapply(1:3, function(k) .tc_pair_ktau(dep$eta[[k]], fam[k]))
    pl <- .tc_pcc_labels(attr(copula, "root"))$short
    dimnames(nat) <- list(NULL, paste0("par", pl))
    dimnames(kt) <- list(NULL, paste0("tau", pl))
    out <- list(loglik = wsum(ll_i), ll_i = ll_i, R = NULL, rho = NULL,
                pair = nat, ktau = kt, rosenblatt = cp$rosenblatt,
                df = cp$df, df_i = dep$df_i, margins = me, U = U, logcop_i = cp$lc)
  } else {
    rho <- cbind(rho12 = dep$r12, rho13 = dep$r13, rho23 = dep$r23,
                 `rho23|1` = dep$p23)
    out <- list(loglik = wsum(ll_i), ll_i = ll_i,
                R = .tc_R(mean(dep$r12), mean(dep$r13), mean(dep$r23)),
                rho = rho, df = cp$df, margins = me, U = U, logcop_i = cp$lc)
    if (length(cp$df) > 1L) {
      # covariate-dependent df: report the sample average, keep the values
      out$df_i <- cp$df
      out$df <- mean(cp$df)
    }
  }
  if (deriv) {
    S <- matrix(0, length(ll_i), idx$npar)
    if (!is.null(mc)) {
      # censoring: the copula term depends on the margins through the lower
      # and upper limits of every coordinate
      for (j in 1:3) {
        dv <- function(nm) {
          me[[j]][[paste0("dl_", nm)]] + cp$du[, j] * mc$dlo[[j]][[nm]] +
            cp$du_up[, j] * mc$dup[[j]][[nm]]
        }
        S[, idx$beta[[j]]] <- X[[j]] * dv("eta")
        S[, idx$extral[[j]]] <- Z[[j]] * dv("extra")
        if (length(idx$shapel[[j]])) S[, idx$shapel[[j]]] <- Q[[j]] * dv("shape")
      }
    } else
    for (j in 1:3) {
      S[, idx$beta[[j]]] <- X[[j]] *
        (me[[j]]$dl_eta + cp$du[, j] * me[[j]]$du_eta)
      S[, idx$extral[[j]]] <- Z[[j]] *
        (me[[j]]$dl_extra + cp$du[, j] * me[[j]]$du_extra)
      if (length(idx$shapel[[j]])) {
        S[, idx$shapel[[j]]] <- Q[[j]] *
          (me[[j]]$dl_shape + cp$du[, j] * me[[j]]$du_shape)
      }
    }
    S[, c(idx$dep, idx$df, idx$ns)] <- cp$score
    if (!is.null(wt)) S <- S * wt
    out$score <- S
  }
  out
}

# C-vine sampler with observation-specific correlations (vectors r12, r13
# and partial correlation p23 of length n).
.tc_rcopula_vec <- function(n, r12, r13, p23, copula, df) {
  E <- matrix(stats::rnorm(3 * n), n, 3)
  Z <- cbind(E[, 1],
             r12 * E[, 1] + sqrt(1 - r12^2) * E[, 2],
             r13 * E[, 1] + sqrt(1 - r13^2) *
               (p23 * E[, 2] + sqrt(1 - p23^2) * E[, 3]))
  if (copula == "gaussian") return(stats::pnorm(Z))
  S <- stats::rchisq(n, df = df)
  stats::pt(Z / sqrt(S / df), df = df)
}

#' Simulate from a trivariate Gaussian or Student-t copula
#'
#' @param n Number of observations.
#' @param R Positive-definite 3 by 3 correlation matrix.
#' @param copula `"gaussian"` or `"student"`.
#' @param df Degrees of freedom of the Student-t copula.
#' @return An `n` by 3 matrix of uniforms.
#' @examples
#' R <- matrix(c(1, .5, .3, .5, 1, .4, .3, .4, 1), 3)
#' U <- tc_rcopula(5, R, "student", df = 5)
#' @export
tc_rcopula <- function(n, R, copula = c("gaussian", "student"), df = 6) {
  copula <- match.arg(copula)
  if (!is.matrix(R) || !all(dim(R) == c(3, 3)) ||
      max(abs(R - t(R))) > 1e-8 || max(abs(diag(R) - 1)) > 1e-8 ||
      min(eigen(R, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
    stop("R must be a positive-definite 3x3 correlation matrix.")
  }
  Z <- matrix(stats::rnorm(3 * n), n, 3) %*% chol(R)
  if (copula == "gaussian") {
    U <- stats::pnorm(Z)
  } else {
    if (!is.numeric(df) || length(df) != 1L || df <= 0) {
      stop("df must be a positive number.")
    }
    S <- stats::rchisq(n, df = df)
    U <- stats::pt(Z / sqrt(S / df), df = df)
  }
  colnames(U) <- c("u1", "u2", "u3")
  U
}

#' Simulate from a trivariate Gaussian or Student-t copula with normal margins
#'
#' Generates copula uniforms and returns standard-normal margins. Intended
#' for reproducible dependence and optimizer smoke tests; use
#' [tc_rcopula()] with quantile functions, or [simulate.tricopula()], for
#' other margins.
#'
#' @inheritParams tc_rcopula
#' @param seed Optional random seed.
#' @return Data frame with `y1`, `y2`, `y3`.
#' @examples
#' head(tc_simulate(10, seed = 1))
#' @export
tc_simulate <- function(n = 500, R = matrix(c(1, .4, .2, .4, 1, .3, .2, .3, 1), 3),
                        copula = c("gaussian", "student"), df = 6, seed = NULL) {
  copula <- match.arg(copula)
  if (!is.null(seed)) set.seed(seed)
  U <- tc_rcopula(n, R, copula, df)
  data.frame(y1 = stats::qnorm(U[, 1]), y2 = stats::qnorm(U[, 2]),
             y3 = stats::qnorm(U[, 3]))
}
