# Pair-copula construction (C-vine rooted at variable 1) -------------------
#
#   c(u1, u2, u3) = c12(u1, u2) * c13(u1, u3) * c23|1(h12(u2 | u1), h13(u3 | u1))
#
# Each pair copula has one dependence predictor eta (unconstrained):
#   gaussian, student : rho = tanh(eta)   (student: df = 2.05 + exp(tau))
#   clayton           : theta = exp(eta) > 0
#   gumbel            : theta = 1 + exp(eta) > 1
# Rotations of Clayton and Gumbel: suffix 90, 180, 270 with
#   c90(u, v) = c(1 - u, v), c180(u, v) = c(1 - u, 1 - v), c270(u, v) = c(u, 1 - v).
# h(v | u) = dC(u, v) / du is the conditional CDF of the second argument
# given the first.

#   frank             : theta = eta (any real number; 0 is independence)
#   joe               : theta = 1 + exp(eta) > 1
# The Frank copula is symmetric in its tails and covers negative dependence
# through the sign of theta, so it has no rotations.
# Vine root (0.10.0): with root r the vine is evaluated on the responses in
# the order (r, a, b), a < b the other two; the pairs are then (r, a),
# (r, b) and (a, b | r), the root being the first argument of the two
# unconditional pair copulas.

.tc_pair_families <- c("gaussian", "student", "clayton", "gumbel",
                       "clayton90", "clayton180", "clayton270",
                       "gumbel90", "gumbel180", "gumbel270",
                       "frank", "joe", "joe90", "joe180", "joe270")

.tc_pair_base <- function(family) sub("(90|180|270)$", "", family)
.tc_pair_rot <- function(family) {
  r <- sub("^[a-z]+", "", family)
  ifelse(nzchar(r), as.integer(r), 0L)
}

# Largest admissible |eta| of a pair copula: 6 on the atanh / log scales,
# 100 for the Frank parameter (identity scale; Kendall's tau 0.96).
.tc_frank_max <- 100
.tc_pair_eta_lim <- function(family) {
  if (.tc_pair_base(family) == "frank") .tc_frank_max else .tc_eta_max
}

# Order of the responses in the vine and labels of the three pairs.
.tc_pcc_order <- function(root = NULL) {
  if (is.null(root)) root <- 1L
  c(as.integer(root), setdiff(1:3, root))
}
.tc_pcc_labels <- function(root = NULL) {
  o <- .tc_pcc_order(root)
  s <- function(i, j) paste0(min(i, j), max(i, j))
  short <- c(s(o[1], o[2]), s(o[1], o[3]), paste0(o[2], o[3], "|", o[1]))
  list(short = short,
       long = paste0("pcc", c(short[1:2], paste0(o[2], o[3], "_given_", o[1]))))
}

# Gauss-Legendre nodes, computed once per session and size.
.tc_gl_cache <- new.env(parent = emptyenv())
.tc_gl <- function(m) {
  key <- as.character(m)
  if (is.null(.tc_gl_cache[[key]])) .tc_gl_cache[[key]] <- .tc_gauss_legendre(m)
  .tc_gl_cache[[key]]
}

# Kendall's tau of the Frank copula, 1 - 4 / theta * (1 - D1(theta)), with
# the Debye function D1(a) = (1 / a) int_0^a t / (exp(t) - 1) dt computed by
# Gauss-Legendre on 8 panels of [0, min(a, 40)] (the integrand is below
# 2e-16 beyond 40). tau(-theta) = -tau(theta).
.tc_frank_tau <- function(theta) {
  ut <- unique(theta)
  gl <- .tc_gl(16L)
  val <- vapply(ut, function(th) {
    a <- abs(th)
    # series near 0 (the closed form cancels there)
    if (a < 0.1) return(th / 9 - th^3 / 900 + th^5 / 52920)
    up <- min(a, 40)
    w <- up / 8
    t <- as.vector(outer(gl$x * w, (0:7) * w, "+"))
    int <- w * sum(rep(gl$w, 8) * t / expm1(t))
    sign(th) * (1 - 4 / a * (1 - int / a))
  }, numeric(1))
  val[match(theta, ut)]
}

# Kendall's tau of the Joe copula:
#   1 + 2 / (2 - theta) * (digamma(2) - digamma(1 + 2 / theta)),
# with limit 1 - trigamma(2) at theta = 2.
.tc_joe_tau <- function(theta) {
  near <- abs(theta - 2) < 1e-8
  th <- ifelse(near, 3, theta)
  ifelse(near, 1 - trigamma(2),
         1 + 2 / (2 - th) * (digamma(2) - digamma(1 + 2 / th)))
}

# Unrotated pair copula: log-density and h-function.
.tc_pair0 <- function(u, v, eta, base, tau = NULL) {
  # numerical domain: |atanh(rho)| <= 6 for the elliptical pairs,
  # theta <= exp(6) + 1 for Clayton, Gumbel and Joe, |theta| <= 100 for Frank
  if (any(abs(eta) > (if (base == "frank") .tc_frank_max else .tc_eta_max))) {
    n <- max(length(u), length(v), length(eta))
    return(list(l = rep(-Inf, n), h = rep(0.5, n)))
  }
  if (base == "gaussian") {
    r <- tanh(eta)
    x <- stats::qnorm(u); y <- stats::qnorm(v)
    l <- -0.5 * log1p(-r^2) - (r^2 * (x^2 + y^2) - 2 * r * x * y) / (2 * (1 - r^2))
    h <- stats::pnorm((y - r * x) / sqrt(1 - r^2))
  } else if (base == "student") {
    r <- tanh(eta)
    nu <- 2.05 + exp(tau)
    x <- stats::qt(u, nu); y <- stats::qt(v, nu)
    l <- lgamma((nu + 2) / 2) + lgamma(nu / 2) - 2 * lgamma((nu + 1) / 2) -
      0.5 * log1p(-r^2) -
      (nu + 2) / 2 * log1p((x^2 + y^2 - 2 * r * x * y) / (nu * (1 - r^2))) +
      (nu + 1) / 2 * (log1p(x^2 / nu) + log1p(y^2 / nu))
    h <- stats::pt((y - r * x) / sqrt((nu + x^2) * (1 - r^2) / (nu + 1)), nu + 1)
  } else if (base == "clayton") {
    th <- exp(eta)
    lu <- -th * log(u); lv <- -th * log(v)
    m <- pmax(lu, lv)
    A <- m + log(exp(lu - m) + exp(lv - m) - exp(-m))  # log(u^-th + v^-th - 1)
    l <- log1p(th) - (1 + th) * (log(u) + log(v)) - (2 + 1 / th) * A
    h <- exp(-(th + 1) * log(u) - (1 / th + 1) * A)
  } else if (base == "gumbel") {
    th <- 1 + exp(eta)
    x <- -log(u); y <- -log(v)
    lx <- log(x); ly <- log(y)
    m <- pmax(th * lx, th * ly)
    lt <- m + log(exp(th * lx - m) + exp(th * ly - m))   # log(x^th + y^th)
    s <- exp(lt / th)                                    # (x^th + y^th)^(1/th)
    l <- -s + x + y + (th - 1) * (lx + ly) + (1 / th - 2) * lt + log(s + th - 1)
    h <- exp(-s + (1 / th - 1) * lt + (th - 1) * lx + x)
  } else if (base == "frank") {
    # c(u, v; -theta) = c(1 - u, v; theta), so only theta > 0 is evaluated.
    # With D = exp(-th u)(1 - exp(-th v)) + exp(-th v)(1 - exp(-th (1 - v)))
    # (two positive terms, no cancellation):
    #   log c = log th + log(1 - exp(-th)) - th (u + v) - 2 log D,
    #   h(v | u) = exp(-th u)(1 - exp(-th v)) / D.
    n <- max(length(u), length(v), length(eta))
    neg <- rep_len(eta < 0, n)
    th <- pmax(abs(eta), 1e-7)
    uu <- ifelse(neg, 1 - u, u)
    t1 <- exp(-th * uu) * (-expm1(-th * v))
    t2 <- exp(-th * v) * (-expm1(-th * (1 - v)))
    D <- t1 + t2
    l <- log(th) + log(-expm1(-th)) - th * (uu + v) - 2 * log(D)
    h <- t1 / D
  } else if (base == "joe") {
    # S = a + b - a b with a = (1 - u)^th, b = (1 - v)^th:
    #   log c = (1 / th - 2) log S + (th - 1) log((1 - u)(1 - v)) + log(th - 1 + S),
    #   h(v | u) = S^(1 / th - 1) (1 - u)^(th - 1) (1 - b).
    th <- 1 + exp(eta)
    l1u <- log1p(-u); l1v <- log1p(-v)
    la <- th * l1u; lb <- th * l1v
    m <- pmax(la, lb)
    lS <- m + log(exp(la - m) + exp(lb - m) - exp(la + lb - m))
    l <- (1 / th - 2) * lS + (th - 1) * (l1u + l1v) + log(th - 1 + exp(lS))
    h <- exp((1 / th - 1) * lS + (th - 1) * l1u) * (-expm1(lb))
  } else {
    stop("Unsupported pair-copula family: ", base)
  }
  list(l = l, h = h)
}

# Pair copula with rotation; h is clipped away from 0 and 1.
.tc_pair <- function(u, v, eta, family, tau = NULL) {
  base <- .tc_pair_base(family)
  rot <- .tc_pair_rot(family)
  uu <- if (rot %in% c(90L, 180L)) 1 - u else u
  vv <- if (rot %in% c(180L, 270L)) 1 - v else v
  r <- .tc_pair0(uu, vv, eta, base, tau)
  h <- if (rot %in% c(180L, 270L)) 1 - r$h else r$h
  list(l = r$l, h = .tc_clip(h))
}

# Pair copula with first derivatives of l and h with respect to u, v, eta
# (and tau for the Student-t pair) by central differences on these scalar
# arguments. Steps in u and v are relative to the distance from the boundary.
.tc_pair_d <- function(u, v, eta, family, tau = NULL) {
  f0 <- .tc_pair(u, v, eta, family, tau)
  cd <- function(fp, fm, step) {
    list(l = (fp$l - fm$l) / (2 * step), h = (fp$h - fm$h) / (2 * step))
  }
  su <- 1e-5 * pmin(u, 1 - u)
  sv <- 1e-5 * pmin(v, 1 - v)
  se <- 1e-5
  du <- cd(.tc_pair(u + su, v, eta, family, tau), .tc_pair(u - su, v, eta, family, tau), su)
  dv <- cd(.tc_pair(u, v + sv, eta, family, tau), .tc_pair(u, v - sv, eta, family, tau), sv)
  de <- cd(.tc_pair(u, v, eta + se, family, tau), .tc_pair(u, v, eta - se, family, tau), se)
  out <- list(l = f0$l, h = f0$h, dl_u = du$l, dh_u = du$h, dl_v = dv$l, dh_v = dv$h,
              dl_eta = de$l, dh_eta = de$h)
  if (!is.null(tau)) {
    dt <- cd(.tc_pair(u, v, eta, family, tau + se), .tc_pair(u, v, eta, family, tau - se), se)
    out$dl_tau <- dt$l; out$dh_tau <- dt$h
  }
  out
}

# Natural parameter and Kendall's tau of a pair copula from its predictor.
.tc_pair_natural <- function(eta, family) {
  switch(.tc_pair_base(family),
         gaussian = tanh(eta), student = tanh(eta),
         clayton = exp(eta), gumbel = 1 + exp(eta),
         frank = eta, joe = 1 + exp(eta))
}
.tc_pair_ktau <- function(eta, family) {
  base <- .tc_pair_base(family)
  sgn <- if (.tc_pair_rot(family) %in% c(90L, 270L)) -1 else 1
  th <- .tc_pair_natural(eta, family)
  sgn * switch(base, gaussian = 2 / pi * asin(th), student = 2 / pi * asin(th),
               clayton = th / (th + 2), gumbel = 1 - 1 / th,
               frank = .tc_frank_tau(th), joe = .tc_joe_tau(th))
}
# Predictor value giving Kendall's tau `kt` (used for starting values).
.tc_pair_eta_from_tau <- function(kt, family) {
  base <- .tc_pair_base(family)
  if (.tc_pair_rot(family) %in% c(90L, 270L)) kt <- -kt
  if (base %in% c("gaussian", "student")) {
    return(atanh(pmin(pmax(sin(pi * kt / 2), -0.99), 0.99)))
  }
  if (base == "frank") {
    kt <- pmin(pmax(kt, -0.9), 0.9)
    return(vapply(kt, function(k) {
      stats::uniroot(function(e) .tc_frank_tau(e) - k, c(-45, 45), tol = 1e-8)$root
    }, numeric(1)))
  }
  kt <- pmin(pmax(kt, 0.03), 0.95)
  if (base == "joe") {
    return(vapply(kt, function(k) {
      stats::uniroot(function(e) .tc_joe_tau(1 + exp(e)) - k, c(-5.5, 5.5),
                     tol = 1e-8)$root
    }, numeric(1)))
  }
  if (base == "clayton") log(2 * kt / (1 - kt)) else log(1 / (1 - kt) - 1)
}
# Lower and upper tail-dependence coefficients of a pair copula.
.tc_pair_taildep <- function(eta, family, tau = NULL) {
  base <- .tc_pair_base(family); rot <- .tc_pair_rot(family)
  th <- .tc_pair_natural(eta, family)
  lu <- switch(base,
    gaussian = c(0, 0),
    student = {
      nu <- 2.05 + exp(tau)
      rep(2 * stats::pt(-sqrt((nu + 1) * (1 - th) / (1 + th)), nu + 1), 2)
    },
    clayton = c(2^(-1 / th), 0),
    gumbel = c(0, 2 - 2^(1 / th)),
    frank = c(0, 0),
    joe = c(0, 2 - 2^(1 / th)))
  if (rot == 180L) lu <- rev(lu)
  if (rot %in% c(90L, 270L)) lu <- c(0, 0)
  stats::setNames(lu, c("lower", "upper"))
}

# Positions of the coefficient blocks inside the dependence vector
# theta = (delta_12, delta_13, delta_23|1, tau for each Student-t pair).
# V: optional list of three design matrices for log(df - 2.05) of the
# Student-t pairs (NULL, or one intercept column: constant df).
.tc_df_const <- function(Vk) {
  is.null(Vk) || (ncol(Vk) == 1L && identical(colnames(Vk), "(Intercept)"))
}
.tc_pcc_layout <- function(W, families, V = NULL) {
  q <- vapply(W, ncol, integer(1))
  ends <- cumsum(q)
  is_t <- .tc_pair_base(families) == "student"
  r <- vapply(1:3, function(k) {
    if (!is_t[k]) 0L else if (is.null(V)) 1L else ncol(V[[k]])
  }, integer(1))
  tend <- ends[3] + cumsum(r)
  list(pos = lapply(1:3, function(k) ends[k] - q[k] + seq_len(q[k])),
       tpos = lapply(1:3, function(k) tend[k] - r[k] + seq_len(r[k])),
       is_t = is_t, npar = tend[3])
}

# C-vine log-density for uniforms U (columns in vine order: root first);
# with deriv = TRUE also the score matrix with respect to theta and
# d logc / d u.
.tc_pcc_eval <- function(theta, U, W, families, deriv = FALSE, V = NULL, ns = NULL) {
  lay <- .tc_pcc_layout(W, families, V)
  eta <- lapply(1:3, function(k) drop(W[[k]] %*% theta[lay$pos[[k]]]))
  if (!is.null(ns)) {
    # non-simplified vine (0.15.0): the predictor of the conditional pair
    # depends on the conditioning value through a spline
    npos <- lay$npar + seq_len(ns$ncoef)
    Bu <- ns$basis(U[, 1], deriv)
    eta[[3]] <- eta[[3]] + drop(Bu$B %*% theta[npos])
  }
  tau <- lapply(1:3, function(k) {
    if (!lay$is_t[k]) return(NULL)
    if (.tc_df_const(V[[k]])) theta[[lay$tpos[[k]]]] else
      as.vector(V[[k]] %*% theta[lay$tpos[[k]]])
  })
  vmul <- function(k, x) if (.tc_df_const(V[[k]])) x else V[[k]] * x
  dep <- list(eta = eta, families = families, tau = tau)
  # df: sample average per pair (the value itself when df is constant)
  df <- vapply(tau, function(t) if (is.null(t)) NA_real_ else mean(2.05 + exp(t)), numeric(1))
  if (any(lengths(tau) > 1L)) {
    dep$df_i <- sapply(tau, function(t) {
      if (is.null(t)) rep(NA_real_, nrow(U)) else rep_len(2.05 + exp(t), nrow(U))
    })
  }
  if (!deriv) {
    p1 <- .tc_pair(U[, 1], U[, 2], eta[[1]], families[1], tau[[1]])
    p2 <- .tc_pair(U[, 1], U[, 3], eta[[2]], families[2], tau[[2]])
    p3 <- .tc_pair(p1$h, p2$h, eta[[3]], families[3], tau[[3]])
    return(list(lc = p1$l + p2$l + p3$l, dep = dep, df = df,
                rosenblatt = cbind(U[, 1], p1$h, p3$h)))
  }
  d1 <- .tc_pair_d(U[, 1], U[, 2], eta[[1]], families[1], tau[[1]])
  d2 <- .tc_pair_d(U[, 1], U[, 3], eta[[2]], families[2], tau[[2]])
  d3 <- .tc_pair_d(d1$h, d2$h, eta[[3]], families[3], tau[[3]])
  S <- matrix(0, nrow(U), length(theta))
  S[, lay$pos[[1]]] <- W[[1]] * (d1$dl_eta + d3$dl_u * d1$dh_eta)
  S[, lay$pos[[2]]] <- W[[2]] * (d2$dl_eta + d3$dl_v * d2$dh_eta)
  S[, lay$pos[[3]]] <- W[[3]] * d3$dl_eta
  if (lay$is_t[1]) S[, lay$tpos[[1]]] <- vmul(1, d1$dl_tau + d3$dl_u * d1$dh_tau)
  if (lay$is_t[2]) S[, lay$tpos[[2]]] <- vmul(2, d2$dl_tau + d3$dl_v * d2$dh_tau)
  if (lay$is_t[3]) S[, lay$tpos[[3]]] <- vmul(3, d3$dl_tau)
  du <- cbind(d1$dl_u + d2$dl_u + d3$dl_u * d1$dh_u + d3$dl_v * d2$dh_u,
              d1$dl_v + d3$dl_u * d1$dh_v,
              d2$dl_v + d3$dl_v * d2$dh_v)
  if (!is.null(ns)) {
    S[, npos] <- Bu$B * d3$dl_eta
    du[, 1] <- du[, 1] + d3$dl_eta * drop(Bu$dB %*% theta[npos])
  }
  list(lc = d1$l + d2$l + d3$l, dep = dep, df = df, score = S, du = du,
       rosenblatt = cbind(U[, 1], d1$h, d3$h))
}

# Starting values for the dependence vector from pairwise Kendall's tau.
.tc_pcc_start <- function(U, W, families, V = NULL) {
  lay <- .tc_pcc_layout(W, families, V)
  # missing responses: starting values from the complete rows
  if (anyNA(U)) U <- U[stats::complete.cases(U), , drop = FALSE]
  kt <- c(stats::cor(U[, 1], U[, 2], method = "kendall"),
          stats::cor(U[, 1], U[, 3], method = "kendall"), NA)
  e1 <- .tc_pair_eta_from_tau(kt[1], families[1])
  e2 <- .tc_pair_eta_from_tau(kt[2], families[2])
  t0 <- log(8 - 2.05)
  a <- .tc_pair(U[, 1], U[, 2], e1, families[1], if (lay$is_t[1]) t0)$h
  b <- .tc_pair(U[, 1], U[, 3], e2, families[2], if (lay$is_t[2]) t0)$h
  kt[3] <- stats::cor(a, b, method = "kendall")
  e <- c(e1, e2, .tc_pair_eta_from_tau(kt[3], families[3]))
  th <- numeric(0)
  for (k in 1:3) {
    v <- numeric(ncol(W[[k]]))
    ic <- match("(Intercept)", colnames(W[[k]]))
    if (!is.na(ic)) v[ic] <- e[k]
    th <- c(th, v)
  }
  if (is.null(V)) return(c(th, rep(t0, sum(lay$is_t))))
  for (k in which(lay$is_t)) {
    v <- numeric(ncol(V[[k]]))
    ic <- match("(Intercept)", colnames(V[[k]]))
    if (!is.na(ic)) v[ic] <- t0
    th <- c(th, v)
  }
  th
}

# Inverse h-function by vectorized bisection: v such that h(v | u) = p.
.tc_pair_hinv <- function(p, u, eta, family, tau = NULL) {
  n <- max(length(p), length(u), length(eta))
  lo <- rep(0, n); hi <- rep(1, n)
  for (i in 1:60) {
    mid <- (lo + hi) / 2
    below <- .tc_pair(u, mid, eta, family, tau)$h < p
    lo[below] <- mid[below]; hi[!below] <- mid[!below]
  }
  .tc_clip((lo + hi) / 2)
}

# C-vine sampler; eta is a list of three predictors (scalars or vectors).
.tc_rpcc <- function(n, eta, families, tau, ns = NULL, gam = NULL) {
  w <- matrix(stats::runif(3 * n), n, 3)
  u1 <- .tc_clip(w[, 1])
  if (!is.null(ns)) eta[[3]] <- eta[[3]] + drop(ns$basis(u1, FALSE)$B %*% gam)
  u2 <- .tc_pair_hinv(w[, 2], u1, eta[[1]], families[1], tau[[1]])
  b <- .tc_pair_hinv(w[, 3], .tc_clip(w[, 2]), eta[[3]], families[3], tau[[3]])
  u3 <- .tc_pair_hinv(b, u1, eta[[2]], families[2], tau[[2]])
  cbind(u1 = u1, u2 = u2, u3 = u3)
}

# Joint CDF C(u1, u2, u3) of the C-vine at one point, using only h-functions:
#   C = int_0^u1 C23|1(a(s), b(s)) ds,  C23|1(a, b) = int_0^a h23(b | t) dt,
# with a(s) = h12(u2 | s), b(s) = h13(u3 | s). Gauss-Legendre in both
# directions.
.tc_pcc_cdf <- function(u, eta, families, tau, gl = .tc_gauss_legendre(48L), ns = NULL,
                        gam = NULL) {
  s <- .tc_clip(u[1] * gl$x)
  e3 <- if (is.null(ns)) eta[3] else {
    rep(eta[3] + drop(ns$basis(s, FALSE)$B %*% gam), times = length(gl$x))
  }
  a <- .tc_pair(s, u[2], eta[1], families[1], tau[[1]])$h
  b <- .tc_pair(s, u[3], eta[2], families[2], tau[[2]])$h
  tt <- .tc_clip(outer(a, gl$x))                  # inner nodes, one row per s
  hh <- .tc_pair(as.vector(tt), rep(b, times = length(gl$x)), e3,
                 families[3], tau[[3]])$h
  inner <- a * drop(matrix(hh, length(s)) %*% gl$w)
  u[1] * sum(gl$w * inner)
}

#' Simulate from a trivariate pair-copula construction
#'
#' C-vine with pair copulas `c12`, `c13` and `c23|1` when the root is the
#' first variable (the default). With `root = 2` the pairs are `c12`
#' (variable 2 as first argument), `c23` and `c13|2`; with `root = 3` they
#' are `c13`, `c23` (variable 3 as first argument) and `c12|3`.
#'
#' @param n Number of observations.
#' @param families Three pair-copula families, in the order root with the
#'   lower-numbered other variable, root with the higher-numbered other
#'   variable, conditional pair; see [tricopula()].
#' @param par Three parameters on the natural scale: a correlation for
#'   `"gaussian"` and `"student"`, `theta > 0` for Clayton, `theta > 1`
#'   for Gumbel and Joe families (also for their rotations) and any nonzero
#'   `theta` for `"frank"`.
#' @param df Degrees of freedom for the Student-t pairs (one value, or one
#'   per pair; ignored for other families).
#' @param root Root variable of the vine: 1, 2 or 3.
#' @return An `n` by 3 matrix of uniforms, columns in the order of the
#'   variables (not of the vine).
#' @examples
#' U <- tc_rpcc(5, c("clayton", "gaussian", "gumbel180"), c(2, 0.4, 1.5))
#' U2 <- tc_rpcc(5, c("frank", "joe", "gaussian"), c(4, 1.8, 0.3), root = 2)
#' @export
tc_rpcc <- function(n, families, par, df = 6, root = 1) {
  if (length(families) != 3L || any(!families %in% .tc_pair_families)) {
    stop("families must be three of: ", paste(.tc_pair_families, collapse = ", "))
  }
  if (!is.numeric(par) || length(par) != 3L) stop("par must be three numbers.")
  if (length(root) != 1L || !root %in% 1:3) stop("root must be 1, 2 or 3.")
  df <- rep_len(df, 3L)
  base <- .tc_pair_base(families)
  eta <- lapply(1:3, function(k) {
    switch(base[k],
      gaussian = , student = {
        if (abs(par[k]) >= 1) stop("correlations must lie in (-1, 1).")
        atanh(par[k])
      },
      clayton = { if (par[k] <= 0) stop("Clayton theta must be positive."); log(par[k]) },
      gumbel = { if (par[k] <= 1) stop("Gumbel theta must exceed 1."); log(par[k] - 1) },
      joe = { if (par[k] <= 1) stop("Joe theta must exceed 1."); log(par[k] - 1) },
      frank = {
        if (par[k] == 0 || abs(par[k]) > .tc_frank_max) {
          stop("Frank theta must be nonzero and at most ", .tc_frank_max, " in absolute value.")
        }
        par[k]
      })
  })
  tau <- lapply(1:3, function(k) if (base[k] == "student") log(df[k] - 2.05) else NULL)
  U <- .tc_rpcc(n, eta, families, tau)
  if (root != 1) {
    U[, .tc_pcc_order(root)] <- U
    colnames(U) <- c("u1", "u2", "u3")
  }
  U
}
