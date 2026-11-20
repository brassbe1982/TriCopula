# Truncation and censoring (0.16.0) -----------------------------------------
#
# Truncation at known bounds [a_j, b_j] of response j, D_j = F_j(b) - F_j(a):
#   type "margin": the margin is the truncated distribution, density
#                  f_j / D_j and PIT (F_j(y) - F_j(a)) / D_j; the copula joins
#                  the truncated margins.
#   type "joint" : the three responses follow the untruncated model and a row
#                  is in the sample only if the truncated response lies inside
#                  its bounds. With one truncated response the density is
#                  c(F_1, F_2, F_3) f_1 f_2 f_3 / D_j: the PIT is not rescaled.
#                  (With two or three truncated responses the divisor is a
#                  rectangle probability of the copula; not implemented.)
#
# Censoring of response j in a row: the response is only known to lie in an
# interval, [y, Inf) (right), (-Inf, y] (left) or [y, y_upper] (interval). On
# the copula scale the censored coordinates K of a row lie in the box
# prod_k [l_k, r_k] and the exact coordinates E are known. The row contributes
#   prod_{j in E} f_j(y_j) * int_box c(u_E, s_K) ds_K,
# the copula density integrated over the censored coordinates. The integral
# is computed by Gauss-Legendre quadrature in every censored coordinate
# (tensor grid for two or three), after the substitution
#   s = l + (r - l) G(x),  G(x) = x^3 (10 - 15 x + 6 x^2),
# whose Jacobian 30 x^2 (1 - x)^2 removes the endpoint singularities of the
# copula density. Only the copula density and its analytic score are needed,
# so every copula of the package is covered, including non-simplified vines.
# Derivatives: with P_q the normalized weights of the nodes,
#   d log I / d theta = sum_q P_q score_q,
#   d log I / d u_e   = sum_q P_q dlogc/du_e (q)              (exact e),
#   d log I / d l_k   = sum_q P_q dlogc/du_k (q) (1 - G_q) - 1 / (r_k - l_k),
#   d log I / d r_k   = sum_q P_q dlogc/du_k (q) G_q + 1 / (r_k - l_k).

# Nodes per censored coordinate for one, two and three censored responses in
# a row (accuracy: tests/testthat/test-censoring.R and the run 15 report).
.tc_cens_nodes <- c(64L, 32L, 16L)

# Check and expand the truncation and censoring arguments of tricopula().
# Returns NULL when neither is used, otherwise a list with
#   lower, upper : per response NULL or vectors of length n (truncation)
#   ttype        : "margin" or "joint"
#   cens         : NULL or n x 3 integer matrix, 0 exact, 1 right, 2 left,
#                  3 interval
#   yup          : per response NULL or the upper limits of interval-censored
#                  values (a valid placeholder elsewhere)
.tc_ct_setup <- function(truncation, truncation_type, censoring, data, y, margins) {
  n <- length(y[[1]])
  lower <- upper <- vector("list", 3)
  if (!is.null(truncation)) {
    if (!is.list(truncation) || length(truncation) != 3L) {
      stop("truncation must be a list of three: NULL, c(lower, upper), or ",
           "list(lower = , upper = ) with numbers or column names of data.")
    }
    bound <- function(b, j, what) {
      if (is.null(b)) return(rep(if (what == "lower") -Inf else Inf, n))
      if (is.character(b) && length(b) == 1L) {
        if (!b %in% names(data)) stop("truncation: column '", b, "' is not in data.")
        b <- data[[b]]
        if (!is.numeric(b) || anyNA(b)) stop("truncation: bounds must be numbers without NA.")
        return(as.numeric(b))
      }
      if (!is.numeric(b) || length(b) != 1L || is.na(b)) {
        stop("truncation: each bound must be one number or one column name.")
      }
      rep(as.numeric(b), n)
    }
    for (j in 1:3) {
      tj <- truncation[[j]]
      if (is.null(tj)) next
      if (is.numeric(tj) && length(tj) == 2L) tj <- list(lower = tj[1], upper = tj[2])
      if (!is.list(tj) || !all(names(tj) %in% c("lower", "upper"))) {
        stop("truncation[[", j, "]] must be NULL, c(lower, upper) or ",
             "list(lower = , upper = ).")
      }
      lo <- bound(tj$lower, j, "lower"); up <- bound(tj$upper, j, "upper")
      if (any(lo >= up)) stop("truncation: lower bounds must be below upper bounds.")
      if (all(is.infinite(lo)) && all(is.infinite(up))) next
      lower[[j]] <- lo; upper[[j]] <- up
    }
  }
  ntr <- sum(!vapply(lower, is.null, logical(1)))
  if (ntr > 1L && truncation_type == "joint") {
    stop("truncation_type = \"joint\" is implemented for one truncated response only.")
  }
  cens <- NULL; yup <- vector("list", 3)
  if (!is.null(censoring)) {
    if (!is.list(censoring) || length(censoring) != 3L) {
      stop("censoring must be a list of three: NULL or c(status = \"column\"",
           "[, upper = \"column\"]).")
    }
    cens <- matrix(0L, n, 3)
    for (j in 1:3) {
      cj <- censoring[[j]]
      if (is.null(cj)) next
      if (!is.character(cj) || is.null(names(cj)) || !"status" %in% names(cj) ||
          !all(names(cj) %in% c("status", "upper")) || !all(cj %in% names(data))) {
        stop("censoring[[", j, "]] must be c(status = \"column\"[, upper = \"column\"]) ",
             "with columns of data.")
      }
      st <- data[[cj[["status"]]]]
      if (anyNA(st) || !all(st %in% 0:3)) {
        stop("censoring: status must be 1 (exact), 0 (right-censored), 2 (left-censored) ",
             "or 3 (interval-censored), without NA.")
      }
      code <- c(1L, 0L, 2L, 3L)[as.integer(st) + 1L]   # 0 exact, 1 right, 2 left, 3 interval
      cens[, j] <- code
      yup[[j]] <- y[[j]]
      if (any(code == 3L)) {
        if (!"upper" %in% names(cj)) {
          stop("censoring[[", j, "]]: interval-censored values need upper = \"column\".")
        }
        up <- data[[cj[["upper"]]]]
        iv <- code == 3L
        if (!is.numeric(up) || anyNA(up[iv]) || any(up[iv] <= y[[j]][iv])) {
          stop("censoring: the upper limit of an interval-censored value must be a ",
               "number above the response.")
        }
        yup[[j]][iv] <- up[iv]
        .tc_check_support(yup[[j]], margins[j])
      }
    }
    if (!any(cens > 0L)) cens <- NULL
  }
  if (!ntr && is.null(cens)) return(NULL)
  for (j in 1:3) {
    if (is.null(lower[[j]])) next
    top <- if (is.null(cens)) y[[j]] else ifelse(cens[, j] == 3L, yup[[j]], y[[j]])
    if (any(y[[j]] < lower[[j]]) || any(top > upper[[j]])) {
      stop("Response ", j, " has values outside its truncation bounds.")
    }
  }
  list(lower = lower, upper = upper, ttype = truncation_type, cens = cens,
       yup = if (is.null(cens)) NULL else yup)
}

# Margin evaluation with truncation: .tc_margin_eval() and the rescaling
# described above. `lo`, `up`: bounds per observation (NULL: no truncation).
.tc_margin_ct <- function(y, X, margin, beta, extra, deriv = FALSE, shape = NULL,
                          lo = NULL, up = NULL, ttype = "margin") {
  m <- .tc_margin_eval(y, X, margin, beta, extra, deriv = deriv, shape = shape)
  if (is.null(lo)) return(m)
  dn <- intersect(c("eta", "extra", "shape"), sub("^du_", "", grep("^du_", names(m), value = TRUE)))
  # CDF and its derivatives at a bound; 0 (or 1) where the bound is outside
  # the support
  sup <- switch(margin, gaussian = c(-Inf, Inf), beta = c(0, 1), c(0, Inf))
  at <- function(b, outside) {
    act <- b > sup[1] & b < sup[2]
    out <- list(F = rep(outside, length(y)))
    out$F[!act & b >= sup[2]] <- 1
    out$F[!act & b <= sup[1]] <- 0
    for (nm in dn) out[[nm]] <- numeric(length(y))
    if (any(act)) {
      yb <- ifelse(act, b, y)
      mb <- .tc_margin_eval(yb, X, margin, beta, extra, deriv = deriv, shape = shape)
      out$F[act] <- mb$u[act]
      for (nm in dn) out[[nm]][act] <- mb[[paste0("du_", nm)]][act]
    }
    out
  }
  a <- at(lo, 0); b <- at(up, 1)
  D <- pmax(b$F - a$F, 1e-300)
  m$logf <- m$logf - log(D)
  # Numerical domain: F(b) - F(a) >= 1e-8. Below that the difference of the
  # two CDF values has lost its digits (and a scale parameter running to
  # infinity makes any truncated margin look flat), so such points are
  # rejected.
  m$logf[D < 1e-8] <- -Inf
  u0 <- m$u
  if (ttype == "margin") m$u <- .tc_clip((u0 - a$F) / D)
  # -log D and its derivatives: the contribution that remains for a censored
  # value under truncation of type "joint"
  m$ld <- -log(D)
  if (deriv) {
    for (nm in dn) {
      dD <- b[[nm]] - a[[nm]]
      m[[paste0("dld_", nm)]] <- -dD / D
      m[[paste0("dl_", nm)]] <- m[[paste0("dl_", nm)]] - dD / D
      if (ttype == "margin") {
        m[[paste0("du_", nm)]] <- (m[[paste0("du_", nm)]] - a[[nm]] - m$u * dD) / D
      }
    }
  }
  m$trunc_prob <- D
  m
}

# The three margins of a model with truncation and/or censoring at parameter
# vector `par`. Returns the margin evaluations `me` (log-density 0 for a
# censored response) and, with censoring, the limits Lo, Up of every
# coordinate on the copula scale, the indicator `cen` and the derivatives of
# the limits (`dlo`, `dup`: per response a list with eta, extra, shape).
.tc_margins_ct <- function(par, y, X, margins, idx, Z, Q, deriv, ct) {
  ev <- function(j, yy) {
    is3 <- idx$shapel[[j]]
    .tc_margin_ct(yy, X[[j]], margins[j], par[idx$beta[[j]]],
                  drop(Z[[j]] %*% par[idx$extral[[j]]]), deriv = deriv,
                  shape = if (length(is3)) drop(Q[[j]] %*% par[is3]),
                  lo = ct$lower[[j]], up = ct$upper[[j]], ttype = ct$ttype)
  }
  me <- lapply(1:3, function(j) ev(j, y[[j]]))
  if (is.null(ct$cens)) return(list(me = me))
  n <- length(y[[1]])
  Lo <- Up <- matrix(0, n, 3)
  dlo <- dup <- vector("list", 3)
  dnm <- c("eta", "extra", "shape")
  for (j in 1:3) {
    code <- ct$cens[, j]
    m <- me[[j]]
    Lo[, j] <- Up[, j] <- m$u
    dj <- intersect(dnm, sub("^du_", "", grep("^du_", names(m), value = TRUE)))
    if (deriv) {
      dlo[[j]] <- lapply(stats::setNames(dj, dj), function(nm) m[[paste0("du_", nm)]])
      dup[[j]] <- lapply(stats::setNames(dj, dj), function(nm) numeric(n))
    }
    if (!any(code > 0L)) next
    rt <- code == 1L; lf <- code == 2L; iv <- code == 3L
    jt <- ct$ttype == "joint" && !is.null(ct$lower[[j]])
    Up[rt, j] <- 1
    Lo[lf, j] <- 0
    if (any(iv)) {
      mu <- ev(j, ct$yup[[j]])
      Up[iv, j] <- mu$u[iv]
    }
    if (deriv) {
      for (nm in dj) {
        d <- m[[paste0("du_", nm)]]
        dlo[[j]][[nm]][lf] <- 0
        dup[[j]][[nm]][lf] <- d[lf]
        if (any(iv)) dup[[j]][[nm]][iv] <- mu[[paste0("du_", nm)]][iv]
        me[[j]][[paste0("dl_", nm)]][code > 0L] <-
          if (jt) m[[paste0("dld_", nm)]][code > 0L] else 0
      }
    }
    me[[j]]$logf[code > 0L] <- if (jt) m$ld[code > 0L] else 0
    me[[j]]$u[code > 0L] <- NA_real_
  }
  list(me = me, Lo = Lo, Up = Up, cen = ct$cens > 0L, dlo = dlo, dup = dup)
}

# Copula contribution with censored coordinates: log of the copula density
# integrated over the box of the censored coordinates (see the top of this
# file). Lo, Up: n x 3 limits (equal for exact coordinates); cen: n x 3
# logical. Value as .tc_cop_eval(), with `du` the derivative with respect to
# the lower limit (the coordinate itself when exact) and `du_up` the
# derivative with respect to the upper limit.
.tc_cop_eval_cens <- function(theta, Lo, Up, cen, W, copula, deriv = FALSE, V = NULL,
                              nodes = .tc_cens_nodes, method = "quadrature") {
  n <- nrow(Lo)
  Um <- Lo
  Um[cen] <- .tc_clip((Lo[cen] + Up[cen]) / 2)
  base <- .tc_cop_eval(theta, Um, W, copula, deriv = FALSE, V = V)
  if (any(!is.finite(base$lc))) return(base)
  lc <- numeric(n)
  S <- matrix(0, n, length(theta))
  dlo <- dup <- matrix(0, n, 3)
  rowsub <- function(M, r) lapply(M, function(m) if (is.null(m)) NULL else m[r, , drop = FALSE])
  pat <- cen[, 1] + 2L * cen[, 2] + 4L * cen[, 3]
  full <- which(pat == 0L)
  if (length(full)) {
    r <- .tc_cop_eval(theta, Lo[full, , drop = FALSE], rowsub(W, full), copula,
                      deriv = deriv, V = if (!is.null(V)) rowsub(V, full))
    lc[full] <- r$lc
    if (deriv) {
      S[full, ] <- r$score
      dlo[full, ] <- r$du
    }
  }
  bad <- function() list(lc = rep(-Inf, n), dep = base$dep, df = base$df)
  for (pt in 1:7) {
    rows <- which(pat == pt)
    if (!length(rows)) next
    Kall <- which(bitwAnd(pt, c(1L, 2L, 4L)) > 0L)
    # method "closed" (0.17.0): one censored coordinate `kin` is integrated
    # in closed form (.tc_cop_inner); quadrature runs over the others, K
    kin <- if (method == "closed") .tc_cens_inner_coord(Kall, copula) else NA_integer_
    K <- if (is.na(kin)) Kall else setdiff(Kall, kin)
    m <- length(K)
    if (m == 0L) {
      r <- .tc_cop_inner(theta, Lo[rows, , drop = FALSE], Lo[rows, kin], Up[rows, kin], kin,
                         rowsub(W, rows), copula, deriv = deriv,
                         V = if (!is.null(V)) rowsub(V, rows))
      if (any(!is.finite(r$lc))) return(bad())
      lc[rows] <- r$lc
      if (deriv) {
        S[rows, ] <- r$score
        dlo[rows, ] <- r$du
        dup[rows, kin] <- r$du_up
      }
      next
    }
    gl <- .tc_gl(nodes[m])
    G <- gl$x^3 * (10 - 15 * gl$x + 6 * gl$x^2)
    gw <- gl$w * 30 * gl$x^2 * (1 - gl$x)^2
    grid <- as.matrix(expand.grid(rep(list(seq_along(G)), m)))
    M <- nrow(grid)
    lgw <- rowSums(matrix(log(gw)[grid], M, m))
    per <- max(1L, floor(50000 / M))
    for (rs in split(rows, ceiling(seq_along(rows) / per))) {
      nr <- length(rs)
      rr <- rep(rs, each = M)
      Ue <- Lo[rr, , drop = FALSE]
      wid <- pmax(Up[rs, K, drop = FALSE] - Lo[rs, K, drop = FALSE], 1e-12)
      for (a in seq_len(m)) {
        Ue[, K[a]] <- .tc_clip(Lo[rr, K[a]] + rep(wid[, a], each = M) * rep(G[grid[, a]], times = nr))
      }
      r <- if (is.na(kin)) {
        .tc_cop_eval(theta, Ue, rowsub(W, rr), copula, deriv = deriv,
                     V = if (!is.null(V)) rowsub(V, rr))
      } else {
        .tc_cop_inner(theta, Ue, Lo[rr, kin], Up[rr, kin], kin, rowsub(W, rr), copula,
                      deriv = deriv, V = if (!is.null(V)) rowsub(V, rr))
      }
      if (any(!is.finite(r$lc))) return(bad())
      lw <- matrix(r$lc, M) + lgw
      mx <- apply(lw, 2, max)
      lse <- mx + log(colSums(exp(sweep(lw, 2, mx))))
      lc[rs] <- lse + rowSums(log(wid))
      if (!deriv) next
      P <- exp(sweep(lw, 2, lse))
      for (j in seq_len(ncol(S))) S[rs, j] <- colSums(P * matrix(r$score[, j], M))
      for (e in setdiff(1:3, K)) dlo[rs, e] <- colSums(P * matrix(r$du[, e], M))
      if (!is.na(kin)) dup[rs, kin] <- colSums(P * matrix(r$du_up, M))
      for (a in seq_len(m)) {
        A <- P * matrix(r$du[, K[a]], M)
        Ga <- G[grid[, a]]
        dlo[rs, K[a]] <- colSums(A * (1 - Ga)) - 1 / wid[, a]
        dup[rs, K[a]] <- colSums(A * Ga) + 1 / wid[, a]
      }
    }
  }
  out <- list(lc = lc, dep = base$dep, df = base$df)
  if (!is.null(base$rosenblatt)) {
    ro <- base$rosenblatt
    ro[pat > 0L, ] <- NA_real_
    out$rosenblatt <- ro
  }
  if (deriv) {
    out$score <- S
    out$du <- dlo
    out$du_up <- dup
  }
  out
}

# Marginal maximum likelihood with truncation and censoring (stage 1 of the
# starting values): exact values contribute the (truncated) log-density,
# censored values log(F(upper) - F(lower)). Started from the fit that
# ignores both.
.tc_fit_margin_ct <- function(j, y, X, margins, Z, Q, w, ct) {
  margin <- margins[j]
  p <- ncol(X[[j]]); e <- ncol(Z[[j]])
  three <- margin %in% .tc_margin_names3
  q3 <- if (three) ncol(Q[[j]]) else 0L
  st <- .tc_fit_margin(y[[j]], X[[j]], margin, Z[[j]], Q[[j]], w)
  # a one-margin index for .tc_margins_ct
  idx1 <- list(beta = list(seq_len(p)), extral = list(p + seq_len(e)),
               shapel = list(p + e + seq_len(q3)))
  ct1 <- list(lower = ct$lower[j], upper = ct$upper[j], ttype = ct$ttype,
              cens = if (is.null(ct$cens)) NULL else ct$cens[, j, drop = FALSE],
              yup = ct$yup[j])
  one <- function(th, deriv) {
    is3 <- idx1$shapel[[1]]
    ev <- function(yy) {
      .tc_margin_ct(yy, X[[j]], margin, th[idx1$beta[[1]]], drop(Z[[j]] %*% th[idx1$extral[[1]]]),
                    deriv = deriv, shape = if (length(is3)) drop(Q[[j]] %*% th[is3]),
                    lo = ct1$lower[[1]], up = ct1$upper[[1]], ttype = "margin")
    }
    m <- ev(y[[j]])
    ll <- m$logf
    dn <- c("eta", "extra", if (three) "shape")
    g <- if (deriv) lapply(stats::setNames(dn, dn), function(nm) m[[paste0("dl_", nm)]])
    if (!is.null(ct1$cens)) {
      code <- ct1$cens[, 1]
      cz <- code > 0L
      if (any(cz)) {
        lo <- ifelse(code == 2L, 0, m$u); up <- ifelse(code == 1L, 1, m$u)
        mu <- NULL
        if (any(code == 3L)) {
          mu <- ev(ct1$yup[[1]])
          up[code == 3L] <- mu$u[code == 3L]
        }
        wd <- pmax(up - lo, 1e-300)
        ll[cz] <- log(wd[cz])
        if (deriv) {
          for (nm in dn) {
            d <- m[[paste0("du_", nm)]]
            dl <- ifelse(code == 2L, 0, d)
            du <- ifelse(code == 1L, 0, d)
            if (!is.null(mu)) du[code == 3L] <- mu[[paste0("du_", nm)]][code == 3L]
            g[[nm]][cz] <- ((du - dl) / wd)[cz]
          }
        }
      }
    }
    list(ll = ll, g = g)
  }
  ws <- function(x) if (is.null(w)) sum(x) else sum(w * x)
  wx <- function(M, x) if (is.null(w)) drop(crossprod(M, x)) else drop(crossprod(M, w * x))
  fn <- function(th) {
    v <- try(suppressWarnings(-ws(one(th, FALSE)$ll)), silent = TRUE)
    if (inherits(v, "try-error") || !is.finite(v)) 1e100 else v
  }
  gr <- function(th) {
    r <- try(suppressWarnings(one(th, TRUE)), silent = TRUE)
    if (inherits(r, "try-error")) return(rep(0, length(th)))
    g <- -c(wx(X[[j]], r$g$eta), wx(Z[[j]], r$g$extra), if (three) wx(Q[[j]], r$g$shape))
    g[!is.finite(g)] <- 0
    g
  }
  opt <- stats::optim(st$par, fn, gr, method = "BFGS",
                      control = list(maxit = 1000, reltol = 1e-10))
  list(par = opt$par, logLik = -opt$value, convergence = opt$convergence)
}

# Lines on truncation and censoring for print methods.
.tc_print_ct <- function(x) {
  if (!is.null(x$truncation)) {
    tx <- vapply(1:3, function(j) {
      tj <- x$truncation[[j]]
      if (is.null(tj)) return(NA_character_)
      if (is.numeric(tj)) tj <- list(lower = tj[1], upper = tj[2])
      paste0("response ", j, " in (", if (is.null(tj$lower)) -Inf else tj$lower, ", ",
             if (is.null(tj$upper)) Inf else tj$upper, ")")
    }, character(1))
    cat("Truncation (type ", x$truncation_type, "): ", paste(tx[!is.na(tx)], collapse = "; "),
        "\n", sep = "")
  }
  if (!is.null(x$censored)) {
    cz <- x$censored
    cat("Censored values per response (right / left / interval): ",
        paste(vapply(1:3, function(j) {
          paste(sum(cz[, j] == 1L), sum(cz[, j] == 2L), sum(cz[, j] == 3L), sep = "/")
        }, character(1)), collapse = ", "),
        "; rows with 1, 2, 3 censored responses: ",
        paste(tabulate(rowSums(cz > 0L), 3L), collapse = ", "), "\n", sep = "")
  }
  invisible(NULL)
}

# Bounds of the truncation for new data: constants are repeated, bounds given
# as columns are taken from `newdata` (the estimation sample when NULL).
.tc_new_trunc <- function(object, newdata, n) {
  tr <- object$truncation
  if (is.null(tr)) return(NULL)
  lapply(1:3, function(j) {
    tj <- tr[[j]]
    if (is.null(tj)) return(NULL)
    get <- function(b, what) {
      if (is.null(b)) return(rep(if (what == "lower") -Inf else Inf, n))
      if (is.character(b)) {
        d <- if (is.null(newdata)) object$data else newdata
        if (!b %in% names(d)) stop("newdata needs the truncation column '", b, "'.")
        return(as.numeric(d[[b]]))
      }
      rep(as.numeric(b), n)
    }
    if (is.numeric(tj)) tj <- list(lower = tj[1], upper = tj[2])
    list(lower = get(tj$lower, "lower"), upper = get(tj$upper, "upper"))
  })
}

# CDF values of the truncation bounds (F(a), F(b)) for new data.
.tc_trunc_F <- function(b, eta, extra, margin, shape) {
  sup <- switch(margin, gaussian = c(-Inf, Inf), beta = c(0, 1), c(0, Inf))
  Fa <- ifelse(b$lower <= sup[1], 0,
               .tc_margin_cdf(pmax(b$lower, sup[1] + 1e-300), eta, extra, margin, shape))
  Fb <- ifelse(b$upper >= sup[2], 1,
               .tc_margin_cdf(pmin(b$upper, sup[2]), eta, extra, margin, shape))
  list(Fa = Fa, Fb = Fb)
}

# Closed forms for one censored coordinate (0.17.0) --------------------------
#
# The copula density integrated over ONE coordinate k between lo and up has a
# closed form:
#   Gaussian:  c_jl(u_j, u_l) * [Phi(a_up) - Phi(a_lo)],
#              a = (z - m) sqrt(i_kk), m = -(i_kj z_j + i_kl z_l) / i_kk, with
#              i the inverse correlation matrix and c_jl the bivariate
#              Gaussian copula density of the other two coordinates;
#   Student-t: c_jl(u_j, u_l) * [T_{df+2}(a_up) - T_{df+2}(a_lo)],
#              a = (z - m) sqrt(i_kk (df + 2) / (df + q_jl)), q_jl the
#              quadratic form of (z_j, z_l);
#   C-vine, k a leaf (root r, other leaf a): with w = h_rk(s | u_r) the
#              substitution removes the pair (r, k):
#              c_ra(u_r, u_a) * [H(w_up) - H(w_lo)], H(w) the conditional
#              distribution function of the pair (a, k | r) given
#              v_a = h_ra(u_a | u_r).
# There is no closed form when the root of a vine is the censored coordinate.
# With two or three censored coordinates one is integrated in closed form and
# quadrature runs over the others (.tc_cop_eval_cens, method "closed").
# Derivatives: with respect to the limits analytically (the copula density at
# the limit over the integral); with respect to the predictors and the exact
# coordinates by central differences of the closed form, as for the pair
# copulas (.tc_pair_d).

# Coordinate to integrate in closed form among the censored ones K (NA: none).
.tc_cens_inner_coord <- function(K, copula) {
  if (copula != "pcc") return(K[length(K)])
  root <- attr(copula, "root")
  lv <- setdiff(K, if (is.null(root)) 1L else root)
  if (length(lv)) lv[length(lv)] else NA_integer_
}

# d C(u, v) / d v = P(U <= u | V = v) of a pair copula. The unrotated families
# are exchangeable, so this is the h-function with the arguments swapped; a
# rotation by 90 degrees becomes one by 270 and the reverse.
.tc_pair_h2 <- function(u, v, eta, family, tau = NULL) {
  rot <- .tc_pair_rot(family)
  fam <- if (rot == 90L) sub("90$", "270", family) else
    if (rot == 270L) sub("270$", "90", family) else family
  .tc_pair(v, u, eta, fam, tau)$h
}

# log(F(b) - F(a)) for a symmetric distribution function, computed from the
# upper tails when a > 0 (no cancellation in the right tail).
.tc_ldiff <- function(a, b, pfun) {
  up <- a > 0
  d <- ifelse(up, pfun(a, lower.tail = FALSE) - pfun(b, lower.tail = FALSE),
              pfun(b) - pfun(a))
  log(pmax(d, 1e-300))
}

# Log of the integral for the Gaussian and Student-t copulas. E: n x 3
# predictors (atanh of rho12, rho13, rho23|1); tau: log(df - 2.05).
.tc_inner_ell <- function(E, tau, U, lo, up, k, copula) {
  if (any(abs(E) > .tc_eta_max)) return(rep(-Inf, nrow(U)))
  d <- .tc_dep_vec(E[, 1], E[, 2], E[, 3])
  r12 <- d$r12; r13 <- d$r13; r23 <- d$r23
  det <- (1 - r12^2) * (1 - r13^2) * (1 - d$p23^2)
  jl <- setdiff(1:3, k); j <- jl[1]; l <- jl[2]
  rjl <- switch(k, r23, r13, r12)
  ikk <- (1 - rjl^2) / det
  rkj <- if (k + j == 3L) r12 else if (k + j == 4L) r13 else r23
  rkl <- if (k + l == 3L) r12 else if (k + l == 4L) r13 else r23
  ikj <- (rkl * rjl - rkj) / det
  ikl <- (rkj * rjl - rkl) / det
  gauss <- copula == "gaussian"
  df <- if (gauss) Inf else 2.05 + exp(tau)
  qf <- if (gauss) stats::qnorm else function(p) stats::qt(p, df)
  zj <- qf(U[, j]); zl <- qf(U[, l])
  m <- -(ikj * zj + ikl * zl) / ikk
  zlo <- qf(lo); zup <- qf(up)
  if (gauss) {
    l2 <- -0.5 * log1p(-rjl^2) - (rjl^2 * (zj^2 + zl^2) - 2 * rjl * zj * zl) / (2 * (1 - rjl^2))
    s <- sqrt(ikk)
    lp <- .tc_ldiff((zlo - m) * s, (zup - m) * s,
                    function(x, lower.tail = TRUE) stats::pnorm(x, lower.tail = lower.tail))
  } else {
    q <- (zj^2 + zl^2 - 2 * rjl * zj * zl) / (1 - rjl^2)
    l2 <- lgamma((df + 2) / 2) + lgamma(df / 2) - 2 * lgamma((df + 1) / 2) -
      0.5 * log1p(-rjl^2) - (df + 2) / 2 * log1p(q / df) +
      (df + 1) / 2 * (log1p(zj^2 / df) + log1p(zl^2 / df))
    s <- sqrt(ikk * (df + 2) / (df + q))
    lp <- .tc_ldiff((zlo - m) * s, (zup - m) * s,
                    function(x, lower.tail = TRUE) stats::pt(x, df + 2, lower.tail = lower.tail))
  }
  l2 + lp
}

# Log of the integral for a C-vine, k a leaf. E: n x 3 predictors of the
# pairs (without the spline of a non-simplified vine, which is added here
# from the root coordinate); taus: list of three (NULL or log(df - 2.05)).
.tc_inner_pcc <- function(E, taus, U, lo, up, k, families, root, ns = NULL, gam = NULL) {
  o <- .tc_pcc_order(root)
  pos <- match(k, o)
  ur <- U[, o[1]]
  e3 <- E[, 3]
  if (!is.null(ns)) e3 <- e3 + drop(ns$basis(ur, FALSE)$B %*% gam)
  lim <- vapply(families, .tc_pair_eta_lim, numeric(1))
  if (any(abs(E[, 1]) > lim[1]) || any(abs(E[, 2]) > lim[2]) || any(abs(e3) > lim[3])) {
    return(rep(-Inf, nrow(U)))
  }
  # values at the limits; a limit of exactly 0 or 1 gives exactly 0 or 1
  at <- function(x, f) {
    h <- f(.tc_clip(x))
    h[x <= 0] <- 0; h[x >= 1] <- 1
    h
  }
  if (pos == 3L) {
    p <- .tc_pair(ur, U[, o[2]], E[, 1], families[1], taus[[1]])
    w <- function(x) at(x, function(s) .tc_pair(ur, s, E[, 2], families[2], taus[[2]])$h)
    H <- function(x) at(x, function(s) .tc_pair(p$h, s, e3, families[3], taus[[3]])$h)
  } else {
    p <- .tc_pair(ur, U[, o[3]], E[, 2], families[2], taus[[2]])
    w <- function(x) at(x, function(s) .tc_pair(ur, s, E[, 1], families[1], taus[[1]])$h)
    H <- function(x) at(x, function(s) .tc_pair_h2(s, p$h, e3, families[3], taus[[3]]))
  }
  p$l + log(pmax(H(w(up)) - H(w(lo)), 1e-300))
}

# Copula density integrated over coordinate k between lo and up (vectors),
# in closed form. Column k of U is ignored. Value as .tc_cop_eval(): lc (log
# of the integral), and with deriv = TRUE the score, `du` (derivative with
# respect to the other coordinates and, in column k, to lo) and `du_up`
# (derivative with respect to up).
.tc_cop_inner <- function(theta, U, lo, up, k, W, copula, deriv = FALSE, V = NULL) {
  n <- nrow(U)
  pcc <- copula == "pcc"
  if (pcc) {
    fam <- attr(copula, "families"); root <- attr(copula, "root"); ns <- attr(copula, "ns")
    lay <- .tc_pcc_layout(W, fam, V)
    pos <- lay$pos
    tau <- lapply(1:3, function(j) {
      if (!lay$is_t[j]) return(NULL)
      if (.tc_df_const(V[[j]])) theta[[lay$tpos[[j]]]] else as.vector(V[[j]] %*% theta[lay$tpos[[j]]])
    })
    npos <- if (!is.null(ns)) lay$npar + seq_len(ns$ncoef)
    gam <- if (!is.null(ns)) theta[npos]
    core <- function(E, tau, U) .tc_inner_pcc(E, tau, U, lo, up, k, fam, root, ns, gam)
  } else {
    q <- vapply(W, ncol, integer(1)); ends <- cumsum(q)
    pos <- lapply(1:3, function(j) ends[j] - q[j] + seq_len(q[j]))
    dfc <- copula != "student" || .tc_df_const(V[[1]])
    dpos <- if (copula == "student") ends[3] + seq_len(if (dfc) 1L else ncol(V[[1]]))
    tau <- if (copula != "student") NULL else if (dfc) theta[[dpos]] else
      as.vector(V[[1]] %*% theta[dpos])
    core <- function(E, tau, U) .tc_inner_ell(E, tau, U, lo, up, k, copula)
  }
  E <- do.call(cbind, lapply(1:3, function(j) drop(W[[j]] %*% theta[pos[[j]]])))
  if (!is.matrix(E)) E <- matrix(E, n, 3)
  l0 <- core(E, tau, U)
  out <- list(lc = l0)
  if (!deriv || any(!is.finite(l0))) return(out)
  h <- 1e-5
  S <- matrix(0, n, length(theta))
  for (j in 1:3) {
    Ep <- E; Em <- E
    Ep[, j] <- Ep[, j] + h; Em[, j] <- Em[, j] - h
    d <- (core(Ep, tau, U) - core(Em, tau, U)) / (2 * h)
    S[, pos[[j]]] <- W[[j]] * d
    if (pcc && j == 3L && !is.null(ns)) {
      S[, npos] <- ns$basis(U[, .tc_pcc_order(root)[1]], FALSE)$B * d
    }
  }
  if (pcc) {
    for (j in which(lay$is_t)) {
      tp <- tau; tm <- tau
      tp[[j]] <- tp[[j]] + h; tm[[j]] <- tm[[j]] - h
      d <- (core(E, tp, U) - core(E, tm, U)) / (2 * h)
      S[, lay$tpos[[j]]] <- if (.tc_df_const(V[[j]])) d else V[[j]] * d
    }
  } else if (copula == "student") {
    d <- (core(E, tau + h, U) - core(E, tau - h, U)) / (2 * h)
    S[, dpos] <- if (dfc) d else V[[1]] * d
  }
  du <- matrix(0, n, 3)
  for (j in setdiff(1:3, k)) {
    st <- 1e-5 * pmin(U[, j], 1 - U[, j])
    Up1 <- U; Um1 <- U
    Up1[, j] <- U[, j] + st; Um1[, j] <- U[, j] - st
    du[, j] <- (core(E, tau, Up1) - core(E, tau, Um1)) / (2 * st)
  }
  # limits: d I / d up = c(u with u_k = up), d I / d lo = -c(u with u_k = lo);
  # no dependence on a limit of exactly 0 or 1
  dens <- function(b) {
    Ub <- U
    Ub[, k] <- .tc_clip(b)
    exp(.tc_cop_eval(theta, Ub, W, copula, deriv = FALSE, V = V)$lc - l0)
  }
  inl <- lo > 0; inu <- up < 1
  if (any(inl)) du[inl, k] <- -dens(lo)[inl]
  dup <- numeric(n)
  if (any(inu)) dup[inu] <- dens(up)[inu]
  out$score <- S; out$du <- du; out$du_up <- dup
  out
}

#' Accuracy of the quadrature of a fit with censored responses
#'
#' Re-evaluates the log-likelihood and its gradient at the estimates with
#' more quadrature nodes and reports what the difference means for the
#' estimates.
#'
#' @details
#' The reference uses the method of the fit (`cens_method`) with `nodes`
#' nodes; the default doubles the nodes of the fit. `shift` is the first-order
#' change of each estimate that the finer quadrature would cause, in units of
#' its standard error: the inverse observed information times the change of
#' the gradient. For a fit with `cens_method = "closed"` the rows are also
#' compared with pure quadrature (column `max_diff_quadrature`), which checks
#' the closed forms themselves; that reference has its own quadrature error.
#' [tricopula()] runs this check after every fit with censored responses
#' (`cens_check = TRUE`), stores the result as `fit$cens_accuracy` and warns
#' if an estimate would move by more than 0.1 standard errors.
#'
#' @param object A `tricopula` fit with censored responses.
#' @param nodes Three node numbers of the reference; default twice those of
#'   the fit (at most 400).
#' @return An object of class `tc_cens_check`: a list with `table` (per
#'   number of censored responses in a row: rows, largest and summed change
#'   of the log-likelihood contributions), `loglik_change`, `shift` (per
#'   parameter, in standard errors; `NULL` without covariance matrix),
#'   `max_shift`, `nodes`, `nodes_reference` and `method`.
#' @export
tc_cens_check <- function(object, nodes = NULL) {
  .tc_check(object)
  if (is.null(object$censored)) stop("The fit has no censored responses.")
  n0 <- object$cens_nodes
  if (is.null(nodes)) nodes <- pmin(2L * n0, 400L)
  if (!is.numeric(nodes) || length(nodes) != 3L || anyNA(nodes) || any(nodes < 4) ||
      any(nodes > 400)) stop("nodes must be three numbers between 4 and 400.")
  meth <- if (is.null(object$cens_method)) "quadrature" else object$cens_method
  ev <- function(nd, method, deriv) {
    y <- object$response
    ct <- attr(y, "ct")
    ct$nodes <- as.integer(nd); ct$method <- method
    attr(y, "ct") <- ct
    .tc_eval(object$coefficients, y, object$design, object$margins, .tc_cop(object),
             object$index, deriv = deriv, W = object$dep_design, Z = object$scale_design,
             V = .tc_V(object), Q = object$shape_design, wt = object$weights,
             obs = object$observed)
  }
  V <- object$vcov
  grad <- !is.null(V) && !anyNA(V)
  e0 <- ev(n0, meth, grad)
  e1 <- ev(nodes, meth, grad)
  eq <- if (meth == "closed") ev(nodes, "quadrature", FALSE)
  ncen <- rowSums(object$censored > 0L)
  w <- if (is.null(object$weights)) rep(1, length(ncen)) else object$weights
  d1 <- e1$ll_i - e0$ll_i
  tab <- do.call(rbind, lapply(1:3, function(m) {
    r <- ncen == m
    data.frame(censored = m, rows = sum(r),
               max_diff = if (any(r)) max(abs(d1[r])) else NA_real_,
               sum_diff = if (any(r)) sum(w[r] * d1[r]) else NA_real_,
               max_diff_quadrature = if (any(r) && !is.null(eq)) max(abs((eq$ll_i - e0$ll_i)[r])) else NA_real_)
  }))
  shift <- NULL
  if (grad) {
    dg <- colSums(e1$score) - colSums(e0$score)
    shift <- stats::setNames(drop(V %*% dg) / sqrt(diag(V)), names(object$coefficients))
  }
  structure(list(table = tab, loglik_change = sum(w * d1), shift = shift,
                 max_shift = if (is.null(shift)) NA_real_ else max(abs(shift)),
                 nodes = n0, nodes_reference = as.integer(nodes), method = meth),
            class = "tc_cens_check")
}

#' @export
print.tc_cens_check <- function(x, ...) {
  cat("Quadrature check (method \"", x$method, "\"): nodes ", paste(x$nodes, collapse = "/"),
      " against ", paste(x$nodes_reference, collapse = "/"), "\n", sep = "")
  tab <- x$table
  if (all(is.na(tab$max_diff_quadrature))) tab$max_diff_quadrature <- NULL
  print(tab, row.names = FALSE, digits = 3)
  cat("Change of the log-likelihood: ", signif(x$loglik_change, 3), "\n", sep = "")
  if (!is.null(x$shift)) {
    cat("Largest first-order change of an estimate:\n  ", signif(x$max_shift, 3),
        " standard errors (", names(x$shift)[which.max(abs(x$shift))], ")\n", sep = "")
  }
  invisible(x)
}
