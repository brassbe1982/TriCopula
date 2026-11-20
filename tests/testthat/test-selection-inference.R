# Inference after family selection (0.14.0): sample splitting and the
# bootstrap with the selection repeated

sim_pcc_sel <- function(n = 300, seed = 201) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n))
  U <- tc_rpcc(n, c("clayton", "gumbel", "frank"), c(1.2, 1.5, 2.5))
  d$y1 <- qnorm(U[, 1], 1 + 0.5 * d$x, 1.2)
  d$y2 <- qlnorm(U[, 2], 0.2 + 0.3 * d$x, 0.5)
  d$y3 <- qgamma(U[, 3], shape = 2.5, scale = exp(0.1 + 0.2 * d$x) / 2.5)
  d
}
fm_sel <- list(y1 ~ x, y2 ~ x, y3 ~ x)

test_that("defaults unchanged: tc_select without split, tc_bootstrap without select", {
  skip_on_cran()
  d <- sim_pcc_sel()
  s <- tc_select(fm_sel, d, hessian = FALSE)
  expect_null(s$split)
  f <- tricopula(fm_sel, d, copula = "pcc", families = s$families, hessian = FALSE)
  expect_identical(coef(s$fit), coef(f))
  expect_identical(s$fit$logLik, f$logLik)
  expect_equal(s$fit$selection$criterion, "AIC")
  expect_null(s$fit$selection$split)
  expect_null(f$selection)
  # the ordinary bootstrap ignores the selection record
  b1 <- tc_bootstrap(s$fit, B = 3, type = "nonparametric", seed = 4)
  b2 <- tc_bootstrap(f, B = 3, type = "nonparametric", seed = 4)
  expect_identical(b1$replicates, b2$replicates)
  expect_null(b1$select)
})

test_that("tc_select(split = ): selection and estimation on disjoint rows", {
  skip_on_cran()
  d <- sim_pcc_sel()
  s <- tc_select(fm_sel, d, split = 0.4, seed = 7)
  expect_equal(length(s$rows_select), 120L)
  expect_equal(sort(c(s$rows_select, s$rows_estimate)), seq_len(nrow(d)))
  expect_equal(s$fit$n, 180L)
  expect_false(is.null(s$fit$vcov))
  # equal to doing the two steps by hand
  s1 <- tc_select(fm_sel, d[s$rows_select, ], hessian = FALSE)
  expect_identical(s$families, s1$families)
  expect_identical(s$table, s1$table)
  expect_identical(coef(s$selection_fit), coef(s1$fit))
  f2 <- tricopula(fm_sel, d[s$rows_estimate, ], copula = "pcc", families = s1$families)
  expect_identical(coef(s$fit), coef(f2))
  expect_identical(vcov(s$fit), vcov(f2))
  # reproducible with the seed, different with another one
  s_again <- tc_select(fm_sel, d, split = 0.4, seed = 7, hessian = FALSE)
  expect_identical(s_again$rows_select, s$rows_select)
  s_other <- suppressWarnings(tc_select(fm_sel, d, split = 0.4, seed = 8, hessian = FALSE))
  expect_false(identical(s_other$rows_select, s$rows_select))
  # weights and the root follow the rows
  w <- rep(c(1, 2), length.out = nrow(d))
  sw <- suppressWarnings(tc_select(fm_sel, d, split = 0.5, seed = 9, weights = w, root = "auto",
                                   hessian = FALSE))
  expect_equal(sw$fit$weights, w[sw$rows_estimate])
  expect_equal(sw$fit$root, sw$root)
  expect_equal(nrow(sw$roots), 3L)
  expect_error(tc_select(fm_sel, d, split = 1), "split must be")
  expect_error(tc_select(fm_sel, d, split = c(0.3, 0.5)), "split must be")
  expect_error(tc_select(fm_sel, d[1:30, ], split = 0.5), "fewer than 20")
})

test_that("tc_bootstrap(select = TRUE) repeats the selection in every sample", {
  skip_on_cran()
  d <- sim_pcc_sel()
  s <- tc_select(fm_sel, d)
  bt <- tc_bootstrap(s$fit, B = 5, type = "nonparametric", select = TRUE, seed = 21)
  expect_s3_class(bt, "tc_bootstrap")
  expect_true(bt$select)
  expect_equal(bt$failures, 0L)
  keep <- c(rownames(bt$coefficients), rownames(bt$natural))
  expect_equal(keep, c(grep("^margin", names(coef(s$fit)), value = TRUE),
                       "margin1_sigma", "margin2_sigma", "margin3_shape",
                       "tau12", "tau13", "tau23|1"))
  expect_equal(dim(bt$replicates), c(5L, length(keep)))
  expect_true(all(bt$coefficients$lower <= bt$coefficients$upper))
  expect_equal(bt$natural$model.se, tc_dependence(s$fit)[rownames(bt$natural), "std.error"])
  # shares per pair add to one; the root was not selected, so it stays 1
  expect_equal(as.numeric(tapply(bt$selection$share, bt$selection$pair, sum)[c("12", "13", "23|1")]),
               c(1, 1, 1))
  expect_equal(unname(bt$roots[["1"]]), 1)
  expect_true(bt$same_selection >= 0 && bt$same_selection <= 1)
  expect_output(print(bt), "Selection repeated")
  # independent reconstruction of the first replicate
  set.seed(21)
  ii <- sample.int(nrow(d), replace = TRUE)
  s1 <- tc_select(fm_sel, d[ii, ], hessian = FALSE)
  v1 <- c(coef(s1$fit), setNames(tc_dependence(s1$fit)$estimate,
                                 rownames(tc_dependence(s1$fit))))[keep]
  expect_equal(unname(bt$replicates[1, ]), unname(v1), tolerance = 1e-10)
  expect_equal(bt$selection$share[bt$selection$pair == "12" &
                                    bt$selection$family == s1$families[1]] >= 0.2, TRUE)
  # the parametric bootstrap and a selected root
  bp <- tc_bootstrap(s$fit, B = 3, select = TRUE, seed = 22)
  expect_equal(nrow(bp$replicates), 3L - bp$failures)
  sa <- tc_select(fm_sel, d, root = "auto", hessian = FALSE)
  ba <- tc_bootstrap(sa$fit, B = 4, type = "nonparametric", select = TRUE, seed = 23)
  expect_equal(sum(ba$roots), 1)
  expect_true(all(ba$natural$replicates <= 4L))
  # refusals
  f <- tricopula(fm_sel, d, copula = "pcc", families = s$families, hessian = FALSE)
  expect_error(tc_bootstrap(f, B = 3, select = TRUE), "tc_select")
  ss <- tc_select(fm_sel, d, split = 0.5, seed = 1, hessian = FALSE)
  expect_error(tc_bootstrap(ss$fit, B = 3, select = TRUE), "split")
  expect_error(tc_bootstrap(s$fit, B = 3, select = NA), "TRUE or FALSE")
})

# 0.14.1: tail dependence, means and joint probabilities after selection;
# selection of margins, covariates and copula type

sim_sel2 <- function(n = 300, seed = 202) {
  d <- sim_pcc_sel(n, seed)
  set.seed(seed + 1)
  d$z <- rnorm(n)
  d
}

test_that("0.14.1 defaults: explicit default arguments give the 0.14.0 selection", {
  d <- sim_pcc_sel()
  s <- tc_select(fm_sel, d, hessian = FALSE)
  s2 <- tc_select(fm_sel, d, margins = list("gaussian", "lognormal", "gamma"), copula = "pcc",
                  covariates = "none", hessian = FALSE)
  expect_identical(coef(s$fit), coef(s2$fit))
  expect_identical(s$table, s2$table)
  expect_equal(s$copula, "pcc")
  expect_equal(s$margins, c("gaussian", "lognormal", "gamma"))
  expect_null(s$copulas); expect_null(s$margin_table); expect_null(s$covariate_table)
  f <- tricopula(fm_sel, d, copula = "pcc", families = s$families, hessian = FALSE)
  expect_error(tc_bootstrap(f, B = 2, thresholds = c(0, 1, 1)), "only used with select")
  expect_error(tc_select(fm_sel, d, copula = "vine"), "copula must be")
  expect_error(tc_select(fm_sel, d, margins = list("gaussian", "normal", "gamma")), "margins must be")
})

test_that("quantities common to all selections agree with independent formulas", {
  skip_if_not_installed("mvtnorm")
  d <- sim_pcc_sel()
  f <- tricopula(fm_sel, d, copula = "pcc", families = c("clayton", "gumbel", "frank"))
  nd <- data.frame(x = c(-1, 0.5))
  thr <- rbind(c(0, 0.8, 0.5), c(1.5, 1.5, 1.2))
  ct <- TriCopula:::.tc_common_table(f, thr, nd)
  dn <- tc_dependence(f)
  th <- dn[c("par12", "par13"), "estimate"]
  expect_equal(ct[c("tau12", "tau13", "tau23|1"), "estimate"],
               dn[c("tau12", "tau13", "tau23|1"), "estimate"])
  expect_equal(ct[c("tau12", "tau13", "tau23|1"), "std.error"],
               dn[c("tau12", "tau13", "tau23|1"), "std.error"], tolerance = 1e-5)
  expect_equal(ct[c("lambdaL12", "lambdaU12", "lambdaL13", "lambdaU13"), "estimate"],
               c(2^(-1 / th[1]), 0, 0, 2 - 2^(1 / th[2])))
  # delta method by hand: Clayton theta = exp(eta), lambda = 2^(-1 / theta)
  se_eta <- sqrt(diag(vcov(f)))[["pcc12_clayton_log_theta"]]
  expect_equal(ct["lambdaL12", "std.error"], 2^(-1 / th[1]) * log(2) / th[1] * se_eta,
               tolerance = 1e-5)
  expect_equal(ct[paste0("mean", 1:3), "estimate"],
               unname(colMeans(predict(f, nd, type = "mean"))))
  expect_equal(ct[c("jointprob1", "jointprob2"), "estimate"],
               c(mean(tc_jointprob(f, thr[1, ], nd)), mean(tc_jointprob(f, thr[2, ], nd))))
  expect_true(all(ct[c("jointprob1", "jointprob2"), "std.error"] > 0))
  # elliptical copulas
  g <- tricopula(fm_sel, d, copula = "student", hessian = FALSE)
  cg <- TriCopula:::.tc_common(coef(g), g)
  R <- g$correlation
  expect_equal(unname(cg[c("tau12", "tau13", "tau23")]),
               2 / pi * asin(c(R[1, 2], R[1, 3], R[2, 3])))
  P <- -cov2cor(solve(R))
  expect_equal(unname(cg[c("tau23|1", "tau13|2", "tau12|3")]),
               2 / pi * asin(c(P[2, 3], P[1, 3], P[1, 2])))
  td <- tc_taildep(g)
  expect_equal(unname(cg[c("lambdaL12", "lambdaU13", "lambdaL23")]),
               c(td[1, 2], td[1, 3], td[2, 3]))
})

test_that("selection of margins, covariates and copula type against references", {
  skip_on_cran()
  d <- sim_sel2()
  fm2 <- list(y1 ~ x + z, y2 ~ x + z, y3 ~ x + z)
  s <- suppressWarnings(tc_select(
    fm2, d, margins = list(c("gaussian", "lognormal"), c("lognormal", "gaussian"), "gamma"),
    covariates = "backward", copula = c("gaussian", "pcc"), hessian = FALSE))
  # marginal log-likelihoods from lm(): Gaussian, and log-normal with its Jacobian
  mt <- s$margin_table
  l1g <- as.numeric(logLik(lm(y1 ~ x + z, d)))
  l2g <- as.numeric(logLik(lm(y2 ~ x + z, d)))
  l2l <- as.numeric(logLik(lm(log(y2) ~ x + z, d))) - sum(log(d$y2))
  expect_equal(mt$logLik[mt$margin == 1 & mt$family == "gaussian"], l1g, tolerance = 1e-6)
  expect_true(is.na(mt$logLik[mt$margin == 1 & mt$family == "lognormal"]))  # y1 has negative values
  expect_equal(mt$logLik[mt$margin == 2], c(l2l, l2g), tolerance = 1e-6)
  expect_equal(mt$AIC[mt$margin == 2], -2 * c(l2l, l2g) + 8, tolerance = 1e-6)
  expect_equal(s$margins, c("gaussian", c("lognormal", "gaussian")[which.max(c(l2l, l2g))], "gamma"))
  # backward elimination by AIC against stats::step() for the two lm margins
  st1 <- step(lm(y1 ~ x + z, d), direction = "backward", trace = 0)
  st2 <- step(lm(log(y2) ~ x + z, d), direction = "backward", trace = 0)
  ct <- s$covariate_table
  expect_equal(ct$term[ct$margin == 1 & ct$selected], attr(terms(st1), "term.labels"))
  expect_equal(ct$term[ct$margin == 2 & ct$selected], attr(terms(st2), "term.labels"))
  expect_equal(attr(terms(s$formulas[[1]]), "term.labels"), attr(terms(st1), "term.labels"))
  # copula type: the table holds the criteria of the joint fits
  g <- tricopula(s$formulas, d, margins = s$margins, copula = "gaussian", hessian = FALSE)
  cp <- s$copulas
  expect_equal(cp$AIC[cp$copula == "gaussian"], AIC(g), tolerance = 1e-8)
  expect_equal(cp$copula[cp$selected], cp$copula[which.min(cp$AIC)])
  expect_equal(s$copula, s$fit$copula)
  expect_equal(s$fit$margins, s$margins)
  # an elliptical choice returns that fit
  s_g <- tc_select(fm_sel, d, copula = "gaussian", hessian = FALSE)
  g0 <- tricopula(fm_sel, d, copula = "gaussian", hessian = FALSE)
  expect_identical(coef(s_g$fit), coef(g0))
  expect_null(s_g$families)
  # split: equal to the two steps by hand
  sp <- suppressWarnings(tc_select(
    fm2, d, margins = list("gaussian", c("lognormal", "gamma"), "gamma"),
    covariates = "backward", split = 0.5, seed = 5, hessian = FALSE))
  s1 <- suppressWarnings(tc_select(
    fm2, d[sp$rows_select, ], margins = list("gaussian", c("lognormal", "gamma"), "gamma"),
    covariates = "backward", hessian = FALSE))
  f2 <- tricopula(s1$formulas, d[sp$rows_estimate, ], margins = s1$margins, copula = "pcc",
                  families = s1$families, hessian = FALSE)
  expect_identical(coef(sp$fit), coef(f2))
  expect_identical(sp$margins, s1$margins)
})

test_that("bootstrap with selection: derived quantities, conditional parameters, all steps", {
  skip_on_cran()
  d <- sim_sel2()
  s <- tc_select(fm_sel, d)
  nd <- data.frame(x = 0)
  thr <- c(1, 1.2, 1)
  bt <- tc_bootstrap(s$fit, B = 4, type = "nonparametric", select = TRUE, seed = 31,
                     thresholds = thr, newdata = nd)
  expect_equal(rownames(bt$derived),
               c("lambdaL12", "lambdaL13", "lambdaU12", "lambdaU13", "mean1", "mean2", "mean3",
                 "jointprob1"))
  expect_equal(bt$derived["jointprob1", "estimate"], tc_jointprob(s$fit, thr, nd))
  expect_equal(bt$derived$model.se,
               TriCopula:::.tc_common_table(s$fit, matrix(thr, 1), nd)[rownames(bt$derived), "std.error"])
  # independent reconstruction of the first replicate
  set.seed(31)
  ii <- sample.int(nrow(d), replace = TRUE)
  s1 <- tc_select(fm_sel, d[ii, ], hessian = FALSE)
  expect_equal(bt$derived_replicates[[1, "jointprob1"]], tc_jointprob(s1$fit, thr, nd),
               tolerance = 1e-10)
  expect_equal(unname(bt$derived_replicates[1, c("mean1", "mean2", "mean3")]),
               unname(unlist(predict(s1$fit, nd, type = "mean"))), tolerance = 1e-10)
  td <- tc_taildep(s1$fit)
  expect_equal(unname(bt$derived_replicates[1, c("lambdaL12", "lambdaL13", "lambdaU12", "lambdaU13")]),
               c(td$lower, td$upper), tolerance = 1e-10)
  # conditional parameters: only replicates with the same family count
  expect_equal(rownames(bt$conditional)[1:3], c("par12", "par13", "par23|1"))
  sh <- bt$selection
  for (k in 1:3) {
    pr <- c("12", "13", "23|1")[k]
    expect_equal(bt$conditional$replicates[k],
                 round(4 * sum(sh$share[sh$pair == pr & sh$original])))
  }
  expect_output(print(bt), "Tail dependence")
  # margins, covariates and copula type repeated
  fm2 <- list(y1 ~ x + z, y2 ~ x + z, y3 ~ x + z)
  mg <- list("gaussian", c("lognormal", "gamma"), "gamma")
  s2 <- suppressWarnings(tc_select(fm2, d, margins = mg, covariates = "backward",
                                   copula = c("gaussian", "pcc")))
  b2 <- tc_bootstrap(s2$fit, B = 4, type = "nonparametric", select = TRUE, seed = 32)
  expect_false(any(grepl("^margin2", c(rownames(b2$coefficients), rownames(b2$natural)))))
  expect_true(all(c("margin1_z", "margin3_z") %in% rownames(b2$coefficients)))
  expect_equal(sum(b2$copulas$share), 1)
  expect_equal(sum(b2$margin_selection$share), 1)
  expect_equal(nrow(b2$covariate_selection), 6L)
  expect_true(all(b2$covariate_selection$share >= 0 & b2$covariate_selection$share <= 1))
  # a term that was not selected is a zero in the replicates
  zshare <- b2$covariate_selection$share[b2$covariate_selection$margin == 1 &
                                           b2$covariate_selection$term == "z"]
  expect_equal(mean(b2$replicates[, "margin1_z"] != 0), zshare)
  set.seed(32)
  ii <- sample.int(nrow(d), replace = TRUE)
  r1 <- suppressWarnings(tc_select(fm2, d[ii, ], margins = mg, covariates = "backward",
                                   copula = c("gaussian", "pcc"), hessian = FALSE))
  expect_equal(unname(b2$replicates[1, "margin1_x"]), unname(coef(r1$fit)["margin1_x"]),
               tolerance = 1e-10)
  expect_equal(unname(b2$derived_replicates[1, "mean2"]),
               mean(predict(r1$fit, d, type = "mean")[[2]]), tolerance = 1e-10)
})
