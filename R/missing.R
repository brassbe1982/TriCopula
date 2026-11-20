# Prior weights and missing responses (0.13.0) -------------------------------
#
# Missing responses are handled by the observed-data likelihood, valid under
# missing at random: a row contributes the densities of its observed
# responses and the copula of those responses.
#   three observed : the trivariate copula density
#   two observed   : the bivariate margin of the copula
#   one observed   : no copula term
# Gaussian and Student-t copulas: the bivariate margin is the same family
# with the corresponding correlation (and the same df).
# PCC, in vine order (root, a, b):
#   b missing   : pair copula c(root, a);   a missing : pair copula c(root, b)
#   root missing: c_ab(u_a, u_b) = int_0^1 c_ra(s, u_a) c_rb(s, u_b)
#                 c_ab|r(h(u_a | s), h(u_b | s)) ds, by Gauss-Legendre in s.
#                 Its derivatives are the averages of the full-vine scores at
#                 the nodes, weighted by the integrand.

# Number of Gauss-Legendre nodes for the integral over a missing vine root.
# The nodes are mapped by s = x^3 (10 - 15 x + 6 x^2), whose Jacobian
# 30 x^2 (1 - x)^2 removes the endpoint singularities of the integrand. With
# an all-Gaussian vine, where the integral is known, the result agrees with
# the bivariate Gaussian copula to 1e-10 (1e-4 without the mapping); for
# strongly dependent vines (Clayton 6 / Gumbel 5, Joe 4 / Clayton 4,
# Student-t 0.9 with 3 df) 128 nodes agree with adaptive integration to
# 2e-8 (64 nodes: 3e-4).
.tc_root_nodes <- 128L

# Data preparation: which rows are used and which responses are observed.
# `vars_y`: variables of the three responses; `vars_x`: all other variables
# of the model. mode "complete" keeps complete rows; "likelihood" keeps rows
# with complete covariates and at least one observed response.
.tc_missing_rows <- function(data, vars_y, vars_x, mode) {
  na_x <- if (length(vars_x)) !stats::complete.cases(data[vars_x]) else rep(FALSE, nrow(data))
  na_y <- vapply(vars_y, function(v) {
    if (length(v)) !stats::complete.cases(data[v]) else rep(FALSE, nrow(data))
  }, logical(nrow(data)))
  na_y <- matrix(na_y, nrow(data), 3)
  if (mode == "complete") !na_x & rowSums(na_y) == 0L else !na_x & rowSums(na_y) < 3L
}

# Copula log-density, scores and d logc / du for uniforms with NA (missing
# responses). Same value structure as .tc_cop_eval; rows with one observed
# response contribute 0, `rosenblatt` is NA for incomplete rows.
.tc_cop_eval_missing <- function(theta, U, W, copula, deriv = FALSE, V = NULL) {
  n <- nrow(U)
  obs <- !is.na(U)
  Uf <- U
  Uf[!obs] <- 0.5
  # per-row dependence quantities (and the numerical-domain check) from a
  # plain evaluation on filled uniforms
  base <- .tc_cop_eval(theta, Uf, W, copula, deriv = FALSE, V = V)
  if (any(!is.finite(base$lc))) return(base)
  pcc <- copula == "pcc"
  ord <- if (pcc) .tc_pcc_order(attr(copula, "root")) else 1:3
  ov <- obs[, ord, drop = FALSE]                 # observed, in vine order
  Uv <- U[, ord, drop = FALSE]
  lc <- numeric(n)
  S <- matrix(0, n, length(theta))
  duv <- matrix(0, n, 3)                         # d logc / du in vine order
  rowsub <- function(M, r) lapply(M, function(m) if (is.null(m)) NULL else m[r, , drop = FALSE])
  full <- which(rowSums(obs) == 3L)
  if (length(full)) {
    r <- .tc_cop_eval(theta, U[full, , drop = FALSE], rowsub(W, full), copula,
                      deriv = deriv, V = if (!is.null(V)) rowsub(V, full))
    lc[full] <- r$lc
    if (deriv) {
      S[full, ] <- r$score
      duv[full, ] <- r$du[, ord, drop = FALSE]
    }
  }
  pat <- list(c(1L, 2L), c(1L, 3L), c(2L, 3L))   # observed pair, vine positions
  if (pcc) {
    fam <- attr(copula, "families")
    lay <- .tc_pcc_layout(W, fam, V)
    eta <- base$dep$eta; tau <- base$dep$tau
    at <- function(x, r) if (is.null(x)) NULL else if (length(x) == 1L) x else x[r]
    vmul <- function(k, r, x) if (.tc_df_const(V[[k]])) x else V[[k]][r, , drop = FALSE] * x
    for (k in 1:2) {                             # root observed: one pair copula
      rows <- which(ov[, pat[[k]][1]] & ov[, pat[[k]][2]] & rowSums(ov) == 2L)
      if (!length(rows)) next
      a <- pat[[k]][2]
      if (!deriv) {
        lc[rows] <- .tc_pair(Uv[rows, 1], Uv[rows, a], at(eta[[k]], rows), fam[k],
                             at(tau[[k]], rows))$l
        next
      }
      d <- .tc_pair_d(Uv[rows, 1], Uv[rows, a], at(eta[[k]], rows), fam[k], at(tau[[k]], rows))
      lc[rows] <- d$l
      S[rows, lay$pos[[k]]] <- W[[k]][rows, , drop = FALSE] * d$dl_eta
      if (lay$is_t[k]) S[rows, lay$tpos[[k]]] <- vmul(k, rows, d$dl_tau)
      duv[rows, 1] <- d$dl_u
      duv[rows, a] <- d$dl_v
    }
    rows <- which(!ov[, 1] & ov[, 2] & ov[, 3])  # root missing: integrate it out
    if (length(rows)) {
      gl <- .tc_gl(.tc_root_nodes)
      m <- length(gl$x)
      rr <- rep(rows, each = m)
      sx <- gl$x^3 * (10 - 15 * gl$x + 6 * gl$x^2)
      gw <- gl$w * 30 * gl$x^2 * (1 - gl$x)^2
      Ue <- cbind(.tc_clip(rep(sx, times = length(rows))), Uv[rr, 2], Uv[rr, 3])
      r <- .tc_pcc_eval(theta, Ue, rowsub(W, rr), fam, deriv = deriv,
                        V = if (!is.null(V)) rowsub(V, rr), ns = attr(copula, "ns"))
      lw <- matrix(r$lc, m) + log(gw)            # one column per row
      mx <- apply(lw, 2, max)
      lcr <- mx + log(colSums(exp(sweep(lw, 2, mx))))
      lc[rows] <- lcr
      if (deriv) {
        P <- exp(sweep(lw, 2, lcr))              # weights of the nodes, sum to 1
        for (j in seq_len(ncol(S))) S[rows, j] <- colSums(P * matrix(r$score[, j], m))
        duv[rows, 2] <- colSums(P * matrix(r$du[, 2], m))
        duv[rows, 3] <- colSums(P * matrix(r$du[, 3], m))
      }
    }
  } else {
    q <- vapply(W, ncol, integer(1))
    ends <- cumsum(q)
    pos <- lapply(1:3, function(k) ends[k] - q[k] + seq_len(q[k]))
    dep <- base$dep
    ex <- function(x) rep_len(x, n)
    r12 <- ex(dep$r12); r13 <- ex(dep$r13); r23 <- ex(dep$r23); p23 <- ex(dep$p23)
    sdep <- ex(dep$s)
    rr <- list(r12, r13, r23)
    fam <- if (copula == "student") "student" else "gaussian"
    tau <- if (copula == "student") ex(log(base$df - 2.05)) else NULL
    dfc <- copula != "student" || .tc_df_const(V[[1]])
    dpos <- if (copula == "student") ends[3] + seq_len(if (dfc) 1L else ncol(V[[1]]))
    for (k in 1:3) {
      a <- pat[[k]][1]; b <- pat[[k]][2]
      rows <- which(obs[, a] & obs[, b] & rowSums(obs) == 2L)
      if (!length(rows)) next
      rk <- rr[[k]][rows]
      tk <- if (is.null(tau)) NULL else tau[rows]
      if (!deriv) {
        lc[rows] <- .tc_pair(U[rows, a], U[rows, b], atanh(rk), fam, tk)$l
        next
      }
      d <- .tc_pair_d(U[rows, a], U[rows, b], atanh(rk), fam, tk)
      lc[rows] <- d$l
      duv[rows, a] <- d$dl_u
      duv[rows, b] <- d$dl_v
      dr <- d$dl_eta / (1 - rk^2)                # d logc / d rho_ab
      if (k == 1L) {
        S[rows, pos[[1]]] <- W[[1]][rows, , drop = FALSE] * (dr * (1 - r12[rows]^2))
      } else if (k == 2L) {
        S[rows, pos[[2]]] <- W[[2]][rows, , drop = FALSE] * (dr * (1 - r13[rows]^2))
      } else {
        # rho23 = rho12 rho13 + s rho23|1 depends on all three predictors
        a12 <- r12[rows]; a13 <- r13[rows]; ap <- p23[rows]
        ds12 <- -a12 * sqrt((1 - a13^2) / (1 - a12^2))
        ds13 <- -a13 * sqrt((1 - a12^2) / (1 - a13^2))
        S[rows, pos[[1]]] <- W[[1]][rows, , drop = FALSE] *
          (dr * (a13 + ap * ds12) * (1 - a12^2))
        S[rows, pos[[2]]] <- W[[2]][rows, , drop = FALSE] *
          (dr * (a12 + ap * ds13) * (1 - a13^2))
        S[rows, pos[[3]]] <- W[[3]][rows, , drop = FALSE] *
          (dr * sdep[rows] * (1 - ap^2))
      }
      if (copula == "student") {
        S[rows, dpos] <- if (dfc) d$dl_tau else V[[1]][rows, , drop = FALSE] * d$dl_tau
      }
    }
  }
  out <- list(lc = lc, dep = base$dep, df = base$df)
  if (pcc) {
    ro <- base$rosenblatt
    ro[rowSums(obs) < 3L, ] <- NA_real_
    out$rosenblatt <- ro
  }
  if (deriv) {
    du <- matrix(0, n, 3)
    du[, ord] <- duv
    du[!obs] <- 0
    out$score <- S
    out$du <- du
  }
  out
}
