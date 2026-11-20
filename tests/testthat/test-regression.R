# Numeric regression tests against results frozen with TriCopula 0.8.0
# (validation/make_fixtures.R in the development workspace). A refactor must
# not change any of these numbers; if a test here fails, the model, the
# likelihood or the default estimation path has changed.

fx <- readRDS(test_path("fixtures", "regression_0.8.0.rds"))

fit_fixture <- function(cs) {
  a <- cs$args
  for (nm in intersect(c("scale", "dependence"), names(a))) {
    a[[nm]] <- stats::as.formula(a[[nm]])
  }
  do.call(tricopula, c(list(list(y1 ~ x, y2 ~ x, y3 ~ x), cs$data), a))
}

for (nm in names(fx$cases)) {
  local({
    case_name <- nm
    cs <- fx$cases[[case_name]]
    e <- cs$expected
    f <- fit_fixture(cs)

    test_that(paste0("0.8.0 likelihood and score are unchanged at a fixed point: ",
                     case_name), {
      expect_identical(names(coef(f)), e$par_names)
      expect_equal(tc_loglik(f, cs$par0), e$logLik_par0, tolerance = 1e-10)
      # PCC scores use finite differences for the pair copulas
      expect_equal(unname(tc_gradient(f, cs$par0)), e$gradient_par0,
                   tolerance = 1e-8)
    })

    test_that(paste0("0.8.0 default fit is reproduced to 1e-6: ", case_name), {
      expect_equal(f$convergence, 0)
      expect_lt(max(abs(unname(coef(f)) - e$coefficients)), 1e-6)
      expect_lt(abs(f$logLik - e$logLik), 1e-6)
      expect_lt(max(abs(unname(f$ifm$coefficients) - e$ifm_coefficients)), 1e-6)
      expect_lt(abs(f$ifm$logLik - e$ifm_logLik), 1e-6)
      # Standard errors come from a numerically differentiated score (finite
      # differences for the pair copulas), accurate to about 3e-4 relative:
      # on one machine a refit from a slightly shifted start moves them by up
      # to 2.6e-4 relative while the coefficients agree to 1e-7, and other
      # platforms differ by the same amount. They are therefore compared with
      # a relative tolerance of 1e-3; estimates stay at 1e-6.
      expect_lt(max(abs(unname(sqrt(diag(vcov(f)))) - e$se) / e$se), 1e-3)
      expect_lt(max(abs(unname(sqrt(diag(vcov(f, type = "sandwich")))) -
                          e$se_sandwich) / e$se_sandwich), 1e-3)
      expect_lt(max(abs(unname(as.matrix(utils::head(predict(f, type = "mean"), 5))) -
                          e$fitted_mean_head)), 1e-6)
      expect_lt(max(abs(tc_dependence(f)$estimate - e$dependence)), 1e-5)
      expect_lt(max(abs(unname(c(tc_tau(f))) - e$tau)), 1e-6)
    })
  })
}
