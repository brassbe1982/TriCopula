# Marginal distributions -------------------------------------------------
#
# Each margin has a linear predictor eta = X beta and a second linear
# predictor "extra" = Z gamma on the log scale:
#   gaussian   y ~ N(eta, sigma^2),                    sigma = exp(extra)
#   lognormal  log(y) ~ N(eta, sigma^2),               sigma = exp(extra)
#   gamma      mean exp(eta), shape exp(extra)
#   weibull    scale exp(eta), shape exp(extra):  F(y) = 1 - exp(-(y/scale)^shape)
#   invgauss   mean mu = exp(eta), sigma = exp(extra):
#              f(y) = (2 pi sigma^2 y^3)^(-1/2) exp(-(y - mu)^2 / (2 mu^2 sigma^2 y))
#   beta       mean mu = plogis(eta), precision phi = exp(extra):
#              y ~ Beta(mu phi, (1 - mu) phi)
# The first three are the 0.1.0 parameterizations and are unchanged.
#
# 0.11.0:
#   gpareto    generalized Pareto for exceedances y > 0: scale sigma = exp(eta),
#              shape xi = exp(extra) - 0.5 > -0.5 (the range in which maximum
#              likelihood is regular):
#              F(y) = 1 - (1 + xi y / sigma)^(-1 / xi), y < -sigma / xi if xi < 0
# Three-parameter margins have a third linear predictor "shape" = Q kappa on
# the log scale. With t = (y / b)^a, b = exp(eta), a = exp(extra):
#   dagum      p = exp(shape):  F(y) = (1 + 1 / t)^(-p)
#   singhmaddala  q = exp(shape):  F(y) = 1 - (1 + t)^(-q)

.tc_margin_names <- c("gaussian", "lognormal", "gamma", "weibull", "invgauss", "beta",
                      "gpareto", "dagum", "singhmaddala")
.tc_margin_names3 <- c("dagum", "singhmaddala")

# Generalized Pareto: A(r, xi) = log(1 + xi r) / xi (-> r as xi -> 0), so that
# F = 1 - exp(-A); NaN outside the support. The inverse is r = expm1(xi A) / xi.
.tc_gp_A <- function(r, xi) {
  n <- max(length(r), length(xi))
  r <- rep_len(r, n); xi <- rep_len(xi, n)
  small <- abs(xi) < 1e-8
  A <- suppressWarnings(log1p(xi * r) / xi)
  A[small] <- r[small]
  A
}
.tc_gp_Ainv <- function(A, xi) {
  n <- max(length(A), length(xi))
  A <- rep_len(A, n); xi <- rep_len(xi, n)
  small <- abs(xi) < 1e-8
  r <- expm1(xi * A) / xi
  r[small] <- A[small]
  r
}

# log(1 + exp(x)) without overflow
.tc_softplus <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))

.tc_family_info <- function(margin) {
  switch(margin,
    gaussian = list(extra_label = "log_scale", natural = "sigma", lower = -Inf, upper = Inf, label = "Gaussian"),
    lognormal = list(extra_label = "log_scale", natural = "sigma", lower = 0, upper = Inf, label = "Log-normal"),
    gamma = list(extra_label = "log_shape", natural = "shape", lower = 0, upper = Inf, label = "Gamma"),
    weibull = list(extra_label = "log_shape", natural = "shape", lower = 0, upper = Inf, label = "Weibull"),
    invgauss = list(extra_label = "log_scale", natural = "sigma", lower = 0, upper = Inf, label = "Inverse Gaussian"),
    beta = list(extra_label = "log_precision", natural = "precision", lower = 0, upper = 1, label = "Beta"),
    gpareto = list(extra_label = "log_shape_plus_0.5", natural = "shape_plus_0.5", lower = 0, upper = Inf, label = "Generalized Pareto"),
    dagum = list(extra_label = "log_a", natural = "a", shape_label = "log_p", natural3 = "p", lower = 0, upper = Inf, label = "Dagum"),
    singhmaddala = list(extra_label = "log_a", natural = "a", shape_label = "log_q", natural3 = "q", lower = 0, upper = Inf, label = "Singh-Maddala"),
    stop("Unsupported margin: ", margin))
}

# Numerical guard for inverse-CDF evaluation; this is not data recoding.
.tc_clip <- function(u, eps = 1e-12) pmin(pmax(u, eps), 1 - eps)

.tc_check_support <- function(y, margin) {
  fi <- .tc_family_info(margin)
  if (fi$lower == 0 && fi$upper == Inf && any(y <= 0)) {
    stop(fi$label, " responses must be strictly positive.", call. = FALSE)
  }
  if (fi$upper == 1 && any(y <= 0 | y >= 1)) {
    stop(fi$label, " responses must lie strictly between 0 and 1.", call. = FALSE)
  }
  invisible(TRUE)
}

.tc_pinvgauss <- function(y, mu, sigma) {
  lam <- 1 / sigma^2
  r <- sqrt(lam / y)
  stats::pnorm(r * (y / mu - 1)) +
    exp(2 * lam / mu + stats::pnorm(-r * (y / mu + 1), log.p = TRUE))
}

# Vectorized bisection (the inverse Gaussian quantile has no closed form).
.tc_qinvgauss <- function(p, mu, sigma) {
  n <- max(length(p), length(mu), length(sigma))
  p <- rep_len(p, n); mu <- rep_len(mu, n); sigma <- rep_len(sigma, n)
  lo <- rep(0, n); hi <- mu
  for (i in 1:200) {
    low <- .tc_pinvgauss(hi, mu, sigma) < p
    if (!any(low)) break
    lo[low] <- hi[low]; hi[low] <- 2 * hi[low]
  }
  for (i in 1:80) {
    mid <- (lo + hi) / 2
    below <- .tc_pinvgauss(mid, mu, sigma) < p
    lo[below] <- mid[below]; hi[!below] <- mid[!below]
  }
  (lo + hi) / 2
}

.tc_margin_cdf <- function(q, eta, extra, margin, shape = NULL) {
  switch(margin,
    gpareto = {
      xi <- exp(extra) - 0.5
      r <- q / exp(pmin(eta, 700))
      u <- -expm1(-.tc_gp_A(r, xi))
      u[is.nan(u)] <- 1                    # beyond the upper end point (xi < 0)
      u
    },
    dagum = exp(-exp(shape) * .tc_softplus(-exp(extra) * (log(q) - pmin(eta, 700)))),
    singhmaddala = -expm1(-exp(shape) * .tc_softplus(exp(extra) * (log(q) - pmin(eta, 700)))),
    gaussian = stats::pnorm(q, mean = eta, sd = exp(extra)),
    lognormal = stats::plnorm(q, meanlog = eta, sdlog = exp(extra)),
    gamma = {
      shape <- exp(extra)
      stats::pgamma(q, shape = shape, scale = exp(pmin(eta, 700)) / shape)
    },
    weibull = stats::pweibull(q, shape = exp(extra), scale = exp(pmin(eta, 700))),
    invgauss = .tc_pinvgauss(q, exp(pmin(eta, 700)), exp(extra)),
    beta = {
      mu <- stats::plogis(eta); phi <- exp(extra)
      stats::pbeta(q, mu * phi, (1 - mu) * phi)
    },
    stop("Unsupported margin: ", margin))
}

.tc_margin_quantile <- function(p, eta, extra, margin, shape = NULL) {
  switch(margin,
    gpareto = {
      xi <- exp(extra) - 0.5
      exp(pmin(eta, 700)) * .tc_gp_Ainv(-log1p(-p), xi)
    },
    # t = 1 / (p^(-1 / P) - 1), y = b t^(1 / a)
    dagum = exp(pmin(eta, 700) - log(expm1(-log(p) / exp(shape))) / exp(extra)),
    # t = (1 - p)^(-1 / q) - 1
    singhmaddala = exp(pmin(eta, 700) + log(expm1(-log1p(-p) / exp(shape))) / exp(extra)),
    gaussian = stats::qnorm(p, mean = eta, sd = exp(extra)),
    lognormal = stats::qlnorm(p, meanlog = eta, sdlog = exp(extra)),
    gamma = {
      shape <- exp(extra)
      stats::qgamma(p, shape = shape, scale = exp(pmin(eta, 700)) / shape)
    },
    weibull = stats::qweibull(p, shape = exp(extra), scale = exp(pmin(eta, 700))),
    invgauss = .tc_qinvgauss(p, exp(pmin(eta, 700)), exp(extra)),
    beta = {
      mu <- stats::plogis(eta); phi <- exp(extra)
      stats::qbeta(p, mu * phi, (1 - mu) * phi)
    },
    stop("Unsupported margin: ", margin))
}

# Conditional mean E(Y | x).
.tc_margin_mean <- function(eta, extra, margin, shape = NULL) {
  switch(margin,
    gpareto = {
      xi <- exp(extra) - 0.5
      ifelse(xi < 1, exp(pmin(eta, 700)) / (1 - xi), Inf)
    },
    # b Gamma(p + 1 / a) Gamma(1 - 1 / a) / Gamma(p), finite for a > 1
    dagum = {
      a <- exp(extra); p <- exp(shape)
      ifelse(a > 1, exp(pmin(eta, 700) + lgamma(p + 1 / a) + lgamma(pmax(1 - 1 / a, 1e-300)) -
                          lgamma(p)), Inf)
    },
    # b Gamma(1 + 1 / a) Gamma(q - 1 / a) / Gamma(q), finite for a q > 1
    singhmaddala = {
      a <- exp(extra); q <- exp(shape)
      ifelse(a * q > 1, exp(pmin(eta, 700) + lgamma(1 + 1 / a) +
                              lgamma(pmax(q - 1 / a, 1e-300)) - lgamma(q)), Inf)
    },
    gaussian = eta,
    lognormal = exp(eta + exp(extra)^2 / 2),
    gamma = exp(pmin(eta, 700)),
    weibull = exp(pmin(eta, 700)) * gamma(1 + exp(-extra)),
    invgauss = exp(pmin(eta, 700)),
    beta = stats::plogis(eta),
    stop("Unsupported margin: ", margin))
}

# Distribution parameters on their natural scale (for predict()).
.tc_margin_pars <- function(eta, extra, margin, shape = NULL) {
  switch(margin,
    gpareto = data.frame(scale = exp(pmin(eta, 700)), shape = exp(extra) - 0.5),
    dagum = data.frame(scale = exp(pmin(eta, 700)), a = exp(extra), p = exp(shape)),
    singhmaddala = data.frame(scale = exp(pmin(eta, 700)), a = exp(extra), q = exp(shape)),
    gaussian = data.frame(location = eta, scale = exp(extra)),
    lognormal = data.frame(meanlog = eta, sdlog = exp(extra)),
    gamma = data.frame(mean = exp(pmin(eta, 700)), shape = exp(extra)),
    weibull = data.frame(scale = exp(pmin(eta, 700)), shape = exp(extra)),
    invgauss = data.frame(mean = exp(pmin(eta, 700)), sigma = exp(extra)),
    beta = data.frame(mean = stats::plogis(eta), precision = exp(extra)),
    stop("Unsupported margin: ", margin))
}

# Log-density, PIT and (optionally) first derivatives with respect to the
# linear predictor eta and the unconstrained extra parameter. `extra` may be
# a scalar or a vector with one value per observation. `shape` is the third
# predictor of the three-parameter margins (derivatives dl_shape, du_shape).
.tc_margin_eval <- function(y, X, margin, beta, extra, deriv = FALSE, shape = NULL) {
  eta <- drop(X %*% beta)
  out <- list()
  # Central differences for CDF derivatives without a convenient closed form.
  fd <- function(which) {
    h <- 1e-5
    if (which == "eta") {
      (.tc_margin_cdf(y, eta + h, extra, margin) -
         .tc_margin_cdf(y, eta - h, extra, margin)) / (2 * h)
    } else {
      (.tc_margin_cdf(y, eta, extra + h, margin) -
         .tc_margin_cdf(y, eta, extra - h, margin)) / (2 * h)
    }
  }
  if (margin == "gaussian" || margin == "lognormal") {
    if (margin == "lognormal") .tc_check_support(y, margin)
    sigma <- exp(extra)
    yy <- if (margin == "lognormal") log(y) else y
    z <- (yy - eta) / sigma
    logf <- stats::dnorm(z, log = TRUE) - extra -
      (if (margin == "lognormal") yy else 0)
    u <- stats::pnorm(z)
    pars <- if (margin == "gaussian") list(mu = eta, sigma = sigma) else
      list(meanlog = eta, sdlog = sigma)
    if (deriv) {
      phi <- stats::dnorm(z)
      out$dl_eta <- z / sigma
      out$dl_extra <- z^2 - 1
      out$du_eta <- -phi / sigma
      out$du_extra <- -z * phi
    }
  } else if (margin == "gamma") {
    .tc_check_support(y, margin)
    shape <- exp(extra)
    mu <- exp(pmin(eta, 700))
    logf <- stats::dgamma(y, shape = shape, scale = mu / shape, log = TRUE)
    u <- stats::pgamma(y, shape = shape, scale = mu / shape)
    pars <- list(mean = mu, shape = shape)
    if (deriv) {
      out$dl_eta <- shape * (y / mu - 1)
      out$dl_extra <- shape * (extra + 1 - log(mu) - digamma(shape) +
                                 log(y) - y / mu)
      out$du_eta <- -y * exp(logf)
      # The derivative of the regularized incomplete gamma function with
      # respect to its shape has no closed form: central difference in the
      # log-shape (error O(h^2), about 1e-10).
      out$du_extra <- fd("extra")
    }
  } else if (margin == "weibull") {
    .tc_check_support(y, margin)
    k <- exp(extra)
    eta <- pmin(eta, 700)
    lt <- k * (log(y) - eta)          # log t,  t = (y / scale)^shape
    t <- exp(lt)
    logf <- extra - eta + (k - 1) * (log(y) - eta) - t
    u <- -expm1(-t)
    pars <- list(scale = exp(eta), shape = k)
    if (deriv) {
      out$dl_eta <- k * (t - 1)
      out$dl_extra <- 1 + (1 - t) * lt
      out$du_eta <- -k * t * exp(-t)
      out$du_extra <- t * lt * exp(-t)
    }
  } else if (margin == "invgauss") {
    .tc_check_support(y, margin)
    mu <- exp(pmin(eta, 700))
    sigma <- exp(extra)
    q <- (y - mu)^2 / (mu^2 * sigma^2 * y)
    logf <- -0.5 * log(2 * pi) - extra - 1.5 * log(y) - 0.5 * q
    u <- .tc_pinvgauss(y, mu, sigma)
    pars <- list(mean = mu, sigma = sigma)
    if (deriv) {
      out$dl_eta <- (y - mu) / (mu^2 * sigma^2)
      out$dl_extra <- q - 1
      out$du_eta <- fd("eta")
      out$du_extra <- fd("extra")
    }
  } else if (margin == "beta") {
    .tc_check_support(y, margin)
    mu <- stats::plogis(eta)
    phi <- exp(extra)
    a <- mu * phi; b <- (1 - mu) * phi
    ly <- log(y); l1y <- log1p(-y)
    logf <- lgamma(phi) - lgamma(a) - lgamma(b) + (a - 1) * ly + (b - 1) * l1y
    u <- stats::pbeta(y, a, b)
    pars <- list(mean = mu, precision = phi)
    if (deriv) {
      da <- digamma(a); db <- digamma(b)
      out$dl_eta <- phi * (ly - l1y - da + db) * mu * (1 - mu)
      out$dl_extra <- phi * (digamma(phi) - mu * da - (1 - mu) * db +
                               mu * ly + (1 - mu) * l1y)
      out$du_eta <- fd("eta")
      out$du_extra <- fd("extra")
    }
  } else if (margin == "gpareto") {
    .tc_check_support(y, margin)
    eta <- pmin(eta, 700)
    e <- exp(extra)
    xi <- e - 0.5
    r <- y / exp(eta)
    z <- 1 + xi * r                        # must be positive
    A <- .tc_gp_A(r, xi)
    logf <- -eta - (1 + xi) * A
    logf[!(z > 0)] <- -Inf
    u <- -expm1(-A)
    pars <- list(scale = exp(eta), shape = xi)
    if (deriv) {
      dA_xi <- (r / z - A) / xi
      small <- rep_len(abs(xi) < 1e-8, length(dA_xi))
      dA_xi[small] <- (-r^2 / 2)[small]
      out$dl_eta <- -1 + (1 + xi) * r / z
      out$dl_extra <- e * (-A - (1 + xi) * dA_xi)
      out$du_eta <- -exp(-A) * r / z
      out$du_extra <- exp(-A) * dA_xi * e
    }
  } else if (margin == "dagum" || margin == "singhmaddala") {
    .tc_check_support(y, margin)
    if (is.null(shape)) stop("The ", margin, " margin needs a shape predictor.")
    eta <- pmin(eta, 700)
    a <- exp(extra)
    k <- exp(shape)                        # p (Dagum) or q (Singh-Maddala)
    lt <- a * (log(y) - eta)               # log t, t = (y / b)^a
    w <- stats::plogis(lt)                 # t / (1 + t)
    sp <- .tc_softplus(lt)                 # log(1 + t)
    sm <- .tc_softplus(-lt)                # log(1 + 1 / t)
    if (margin == "dagum") {
      # p log t - (p + 1) log(1 + t) written as -p log(1 + 1/t) - log(1 + t):
      # the first form cancels catastrophically for large p
      logf <- extra + shape - log(y) - k * sm - sp
      u <- exp(-k * sm)
      pars <- list(scale = exp(eta), a = a, p = k)
      if (deriv) {
        dlt <- k * (1 - w) - w             # = p - (p + 1) w without cancellation
        out$dl_eta <- -a * dlt
        out$dl_extra <- 1 + lt * dlt
        out$dl_shape <- 1 - k * sm
        out$du_eta <- -a * u * k * (1 - w)
        out$du_extra <- lt * u * k * (1 - w)
        out$du_shape <- -u * k * sm
      }
    } else {
      logf <- extra + shape - log(y) + lt - (k + 1) * sp
      sv <- exp(-k * sp)                   # survival function
      u <- -expm1(-k * sp)
      pars <- list(scale = exp(eta), a = a, q = k)
      if (deriv) {
        dlt <- 1 - (k + 1) * w
        out$dl_eta <- -a * dlt
        out$dl_extra <- 1 + lt * dlt
        out$dl_shape <- 1 - k * sp
        out$du_eta <- -a * k * w * sv
        out$du_extra <- lt * k * w * sv
        out$du_shape <- k * sp * sv
      }
    }
  } else {
    stop("Unsupported margin: ", margin)
  }
  out$logf <- logf
  out$u <- .tc_clip(u)
  out$pars <- pars
  out
}

# Crude starting values: least squares on a transformed response and a
# moment-type value for the extra parameter.
.tc_margin_start <- function(y, X, margin) {
  yy <- switch(margin, gaussian = y, beta = stats::qlogis(y), log(y))
  fit <- try(stats::lm.fit(X, yy), silent = TRUE)
  bj <- if (inherits(fit, "try-error")) rep(0, ncol(X)) else fit$coefficients
  bj[!is.finite(bj)] <- 0
  s <- max(stats::sd(yy - drop(X %*% bj)), 0.1)
  ic <- match("(Intercept)", colnames(X))
  extra0 <- switch(margin,
    gaussian = log(s), lognormal = log(s),
    gamma = log(2),
    weibull = {
      k <- pi / (sqrt(6) * s)                 # sd of log Weibull = pi / (sqrt(6) k)
      if (!is.na(ic)) bj[ic] <- bj[ic] + 0.5772 / k
      log(k)
    },
    invgauss = log(s / sqrt(mean(y))),         # CV^2 = mu sigma^2
    beta = {
      m <- mean(y)
      log(max(m * (1 - m) / stats::var(y) - 1, 1))
    },
    gpareto = {
      # start at xi = 0.1; near the exponential E log y = log sigma - 0.5772
      if (!is.na(ic)) bj[ic] <- bj[ic] + 0.5772
      log(0.6)
    },
    # p = q = 1 is the log-logistic: sd of log y = pi / (sqrt(3) a)
    dagum = , singhmaddala = log(pi / (sqrt(3) * s)))
  list(beta = bj, extra0 = extra0, shape0 = 0)
}

# Stand-alone marginal maximum likelihood (stage 1 of the IFM start values).
# w: optional weights (prior weight times 0/1 for an observed response).
.tc_fit_margin <- function(y, X, margin, Z, Q = NULL, w = NULL) {
  p <- ncol(X)
  e <- ncol(Z)
  ie <- p + seq_len(e)
  three <- margin %in% .tc_margin_names3
  is3 <- p + e + seq_len(if (three) ncol(Q) else 0L)
  st <- if (is.null(w)) .tc_margin_start(y, X, margin) else
    .tc_margin_start(y[w > 0], X[w > 0, , drop = FALSE], margin)
  ws <- function(x) if (is.null(w)) sum(x) else sum(w * x)
  wx <- function(M, x) if (is.null(w)) drop(crossprod(M, x)) else drop(crossprod(M, w * x))
  # constant scale/shape on the intercept of Z, other coefficients zero
  g0 <- numeric(e)
  ic <- match("(Intercept)", colnames(Z))
  if (!is.na(ic)) g0[ic] <- st$extra0
  k0 <- numeric(length(is3))               # third predictor starts at 0
  sh <- function(th) if (three) drop(Q %*% th[is3]) else NULL
  fn <- function(th) {
    m <- try(suppressWarnings(.tc_margin_eval(y, X, margin, th[seq_len(p)],
                                              drop(Z %*% th[ie]), shape = sh(th))),
             silent = TRUE)
    if (inherits(m, "try-error")) return(1e100)
    v <- -ws(m$logf)
    if (is.finite(v)) v else 1e100
  }
  gr <- function(th) {
    m <- try(suppressWarnings(.tc_margin_eval(y, X, margin, th[seq_len(p)],
                                              drop(Z %*% th[ie]), deriv = TRUE,
                                              shape = sh(th))),
             silent = TRUE)
    if (inherits(m, "try-error")) return(rep(0, length(th)))
    g <- -c(wx(X, m$dl_eta), wx(Z, m$dl_extra), if (three) wx(Q, m$dl_shape))
    g[!is.finite(g)] <- 0
    g
  }
  opt <- stats::optim(c(st$beta, g0, k0), fn, gr, method = "BFGS",
                      control = list(maxit = 1000, reltol = 1e-10))
  list(par = opt$par, logLik = -opt$value, convergence = opt$convergence)
}
