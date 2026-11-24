# TriCopula 0.17.6

Documentation only; no code or test changed.
- The package help page (`?TriCopula`), `README.md` and the vignette list
  the author's applied studies of jointly modelled outcomes (2019 to 2025)
  under "Related publications by the author", with DOIs checked against
  Crossref.

# TriCopula 0.17.5

- Corrected maintainer e-mail address in `DESCRIPTION`:
  Ibrahim.niankara@aau.ac.ae. Nothing else changed.

# TriCopula 0.17.4

No change to the package code: estimates, likelihoods and defaults are
those of 0.17.1.
- `URL` and `BugReports` in `DESCRIPTION` point to the public repository
  <https://github.com/brassbe1982/TriCopula>.
- `README.md`: installation from GitHub.
- Tests: standard errors and covariance matrices are compared with the
  stored results with a relative tolerance of 1e-3 (before: 1e-6 absolute,
  or 1e-6 to 1e-7 relative). The first checks on Linux, macOS and Windows
  Server (GitHub Actions) failed these comparisons only, by at most 2.3e-4,
  while every estimate and log-likelihood agreed to 1e-6 or better. The
  standard errors come from a numerically differentiated score and are
  accurate to about 3e-4 relative: on one machine a refit from a slightly
  shifted start moves them by up to 2.6e-4 relative. The earlier tolerance
  held only because the arithmetic was identical. Estimates and
  likelihoods are still compared at 1e-6 or tighter.
- GitHub Actions: `actions/checkout@v5`.

# TriCopula 0.17.3

Preparation for CRAN; no change to any estimate, likelihood or default.

## Changed
- `DESCRIPTION`: title and description rewritten for CRAN (no
  "prototype" in the title; methods cited with their DOIs); `Depends:
  R (>= 3.5.0)` stated (the test fixtures are serialized in format 3,
  which `R CMD build` otherwise adds silently).
- `print()` of the quadrature check (`fit$cens_accuracy`,
  `tc_cens_check()`) puts the largest first-order change of an estimate on
  a second line, so that the output fits a 70-character console. The
  numbers printed are the same.
- Tests: the 47 slowest of the 107 tests (3.5 seconds or more each, about
  16 minutes together) call `skip_on_cran()`. They run whenever the
  environment variable `NOT_CRAN` is `"true"` (for example
  `Sys.setenv(NOT_CRAN = "true")` before `testthat::test_dir()`, or with
  `devtools::test()`); on CRAN the remaining 60 tests take about one
  minute. No test was removed or weakened. Two tests that compare with
  mvtnorm now call `skip_if_not_installed("mvtnorm")`, so that the tests
  pass when the suggested packages are not installed (checked with a
  library containing only the hard dependencies and testthat).
- The vignette restores the options it sets.

# TriCopula 0.17.2

## Documentation only
- Acknowledgment of the funder, BRASS Digital Lab, in `README.md` and as
  funder (`fnd`) in `DESCRIPTION`.
- No code changed: every result is that of 0.17.1.

# TriCopula 0.17.1

## Changed default
- `cens_method = "closed"` is the default for fits with censored
  responses (it was `"quadrature"` in 0.17.0, the only method of 0.16.0).
  Decided by the package author on the evidence of the run 18 report: fits
  1.4 to 5.7 times faster, quadrature error smaller by two or more orders
  of magnitude, and no fit flagged by the accuracy check under strong
  dependence where a quarter of the quadrature fits were.
  `cens_method = "quadrature"` gives the likelihood of 0.16.0 (tested
  against values stored from 0.16.0). The two methods differ by the
  quadrature error: at most 1e-5 standard errors in the estimates of the
  Gaussian scenario of the Monte Carlo, 0.05 in the vine scenario, more
  under very strong dependence.

## Unchanged
- Everything else; fits without censored responses are not affected.

# TriCopula 0.17.0

Open points of runs 14, 15 and 16 closed (extension run 18).

## New
- `tricopula(cens_method = "closed")`: one censored response per row is
  integrated in closed form and quadrature runs over the remaining ones. A
  row with one censored response then needs no quadrature, a row with two
  needs `cens_nodes[1]` nodes (64 instead of 1024) and a row with three
  `cens_nodes[2]^2` (1024 instead of 4096). Closed forms: conditional
  normal and Student-t distribution functions for the Gaussian and
  Student-t copulas; for a vine the conditional distribution function of
  the conditional pair, for every response except the root (a row in which
  only the root is censored is computed by quadrature as before).
  Non-simplified vines, covariates in every predictor and weights are
  covered. Derivatives with respect to the limits are analytic, the others
  central differences of the closed form.
- `tc_cens_check()` and `tricopula(cens_check = TRUE)`: after a fit with
  censored responses the log-likelihood and its gradient are evaluated
  once more with twice the nodes. The result (`fit$cens_accuracy`) gives
  the change per row and the first-order change of every estimate in
  standard errors; `tricopula()` warns if an estimate would move by more
  than 0.1 standard errors. The object gains `cens_method` and
  `cens_accuracy`.

## Changed default
- `k_cond`, the basis size of the spline of a non-simplified vine, is 10
  instead of 6. Reason (Monte Carlo, n = 1000, 200 replications, cosine
  curve): the lowest pointwise coverage of the 95% interval of the curve
  rises from 0.675 to 0.93, the average from 0.898 to 0.952, and the root
  mean squared error of the curve falls from 0.076 to 0.068; nothing
  changes where the truth is flat or linear. `k_cond = 6` gives the fits
  of 0.15.0 and 0.16.0 (tested against values stored from 0.16.0). A fit
  takes about a quarter longer. Simplified fits are not affected.

## Unchanged
- Every other default. `cens_method = "quadrature"` is the default and
  gives the likelihood of 0.16.0 (tested against values stored from
  0.16.0); `cens_check` adds a component and possibly a warning, not a
  different estimate.

## Limitations
- `cens_method = "closed"` is not the default in this version (it is
  since 0.17.1).
- In the closed forms a difference of two conditional distribution
  functions is formed. For the Gaussian and Student-t copulas it is
  computed from the upper tails where needed; for a vine it loses relative
  accuracy when the conditional probability of the censoring interval is
  below about 1e-8.
- Rows with three censored responses cost about as much as before with
  `cens_method = "closed"` (1024 evaluations with numerical derivatives
  against 4096 with analytic ones) but are more accurate; use smaller
  `cens_nodes` for speed.
- `tc_cens_check()` compares with twice the nodes of the same method. It
  cannot detect an error that both node numbers share.
- The check costs one evaluation of the likelihood and its gradient with
  twice the nodes per coordinate (eight times the work for rows with three
  censored responses under quadrature); `cens_check = FALSE` skips it. The
  bootstrap refits skip it.

# TriCopula 0.16.0

Truncation and censoring (extension run 15).

## New
- `tricopula(truncation = )`: known lower and upper bounds per response,
  as numbers or as columns of the data. `truncation_type = "margin"` (the
  default) makes each truncated response's margin the truncated
  distribution (density `f / (F(b) - F(a))`, rescaled probability integral
  transform); the model describes the truncated population.
  `truncation_type = "joint"` is the sampling case: the untruncated model
  holds and a row is observed only if the truncated response is inside its
  bounds (likelihood divided by `F(b) - F(a)`, transform not rescaled);
  implemented for one truncated response.
- `tricopula(censoring = )`: right, left and interval censoring per
  response, with a status column coded as in `survival::Surv()` (1 exact,
  0 right, 2 left, 3 interval with an `upper` column). A row contributes
  the densities of its exact responses and the copula density integrated
  over its censored ones, by Gauss-Legendre quadrature with analytic
  derivatives (`cens_nodes`, default 64, 32 and 16 nodes per censored
  response for one, two and three censored responses in a row). All
  copulas, vine roots and non-simplified vines; covariates in every
  predictor; weights.
- The starting values use marginal maximum likelihood with the same
  truncation and censoring.
- `predict(type = "mean")`, `simulate()` and `tc_jointprob()` respect the
  truncation; `print()` and `summary()` report truncation and the number
  of censored values. The object gains `truncation`, `truncation_type`,
  `censoring`, `censored`, `cens_nodes`.

## Unchanged
- With `truncation = NULL` and `censoring = NULL` (the defaults) results
  are those of 0.15.0. A status column that is 1 everywhere and bounds
  that are infinite also give identical results.

## Limitations
- Rows with three censored responses cost 4096 evaluations of the copula
  density each, rows with two 1024; fits with many such rows are slow,
  most of all with the Student-t copula.
- The quadrature is accurate to about 1e-6 on the log-likelihood
  contribution at moderate dependence; with Kendall's tau around 0.7 the
  worst rows are off by up to 0.005 (two censored) or 0.2 (three censored)
  with the default nodes. Increase `cens_nodes` to check.
- `truncation_type = "joint"` with two or three truncated responses is
  not implemented.
- Censoring cannot be combined with missing responses. `tc_select()`, the
  bootstrap with selection and the parametric bootstrap (with censoring)
  are not available; the nonparametric bootstrap is.
- `fit$uniforms` and `fit$rosenblatt` are `NA` for censored values, so the
  PIT diagnostics of `tc_diagnostics()` use the exact values only and are
  not valid tests under censoring.
- With `truncation_type = "joint"`, `predict()` gives the mean of the
  untruncated population; `simulate()` draws from the truncated one.
- New suggested package: `survival` (tests only).

# TriCopula 0.15.0

Non-simplified vine (extension run 14).

## New
- `tricopula(copula = "pcc", simplified = FALSE)`: the predictor of the
  conditional pair copula gains a penalized cubic B-spline `f(u_r)` in the
  probability integral transform of the root (`k_cond` basis functions on
  fixed, equally spaced knots, second-order difference penalty, `f`
  integrating to zero). Everything is estimated jointly: `u_r` moves with
  the parameters of the root's margin, and the score contains the
  corresponding chain-rule term. The smoothing parameter is estimated by
  the Fellner-Schall iteration used for `s()` terms, or fixed through `sp`.
  Works with every pair-copula family, every root, covariates in
  `dependence`, weights and missing responses.
- `tc_simplifying_test()`: likelihood-ratio test of the simplified vine
  against an unpenalized spline alternative.
- The spline is listed in `summary()`, returned by `tc_smooth()` (as a
  function of `u`), and used by `simulate()`, `tc_jointprob()`,
  `tc_bootstrap()` and `predict()`.
- New coefficients are the last elements of the parameter vector, named
  e.g. `pcc23_given_1_s(u1).1`; all other names and positions are
  unchanged. The object gains `simplified` and `nonsimplified`.

## Unchanged
- `simplified = TRUE` (the default): results as in 0.14.1.

## Limitations
- Only the conditional pair depends on the conditioning value, and only
  through its dependence parameter (not its degrees of freedom or family).
- `tc_dependence()`, `tc_tau()` and `fit$pair` average over the fitted
  conditioning values; their standard errors treat these values as fixed.
- `predict(type = "dependence")` with `newdata` gives the conditional pair
  at `f = 0` (the average over the conditioning value on the predictor
  scale).
- `tc_select()` selects families for the simplified vine only;
  `tc_bootstrap(select = TRUE)` is not available for non-simplified fits.
- New dependency: package `splines` (part of R).

# TriCopula 0.14.1

Addendum to extension run 16: the cases that 0.14.0 left out.

## New
- `tc_bootstrap(select = TRUE)` also reports `derived`: lower and upper
  tail dependence of the unconditional pairs, the conditional means
  averaged over `newdata`, and joint probabilities for `thresholds` (new
  arguments `thresholds`, `newdata`). These are defined under every
  selection. `conditional` holds the copula parameters from the samples
  that chose the same family (or copula type) as the original fit; these
  are intervals given that choice, not intervals after selection.
- `tc_select()` can also select the family of each margin (`margins` as a
  list of three candidate vectors; marginal fits), the terms of the mean
  formulas (`covariates = "backward"`; marginal fits) and the copula type
  (`copula = "auto"` or a subset of `"gaussian"`, `"student"`, `"pcc"`;
  joint fits). Sample splitting and the bootstrap repeat all steps. The
  result gains `copula`, `margins`, `formulas`, `copulas`, `margin_table`,
  `covariate_table`; the bootstrap gains `copulas`, `margin_selection`,
  `covariate_selection`.

## Unchanged
- With the default arguments (`copula = "pcc"`, `covariates = "none"`,
  three margin names) `tc_select()` and `tc_bootstrap()` give the results
  of 0.14.0; the existing elements of their results are unchanged.

## Limitations
- Coefficients of a margin whose family is selected are not reported by
  the bootstrap (use the means). A term that was not selected counts as a
  coefficient of zero.
- Terms of the `scale`, `shape` and `dependence` formulas are not selected.
  Margin and covariate selection are not available with smooth terms.
- Joint probabilities cost one trivariate integral per distinct row of
  `newdata` and bootstrap sample.

# TriCopula 0.14.0

Inference after selection of the pair-copula families (extension run 16).

## New
- `tc_select(split = , seed = )`: sample splitting. The families (and the
  root, with `root = "auto"`) are selected on a random share `split` of the
  rows; the selected model is estimated on the other rows only. The
  standard errors of the returned fit are then valid given the selection.
  The result also holds `rows_select`, `rows_estimate` and `selection_fit`.
- `tc_bootstrap(select = TRUE)`: the selection is repeated in every
  bootstrap sample (parametric or nonparametric). Reported are the
  quantities that mean the same under every family: coefficients and
  natural parameters of the margins and Kendall's tau of the pairs, with
  bootstrap standard errors and percentile intervals; also `selection`
  (frequency of each family per pair), `roots` and `same_selection`.
  Requires the fit returned by `tc_select()`.
- Fits returned by `tc_select()` carry `fit$selection` (candidates,
  criterion, root setting, further arguments).

## Changed
- `tc_select()`: the preliminary fit that supplies the marginal estimates
  no longer iterates the joint optimizer, whose result was never used. The
  selected families, the selection table and the final fit are identical
  (checked on four models, differences exactly 0).

## Unchanged
- `tc_select()` without `split` and `tc_bootstrap()` without `select`:
  results as in 0.13.1.

## Limitations
- `tc_bootstrap(select = TRUE)` gives no interval for the pair-copula
  parameters themselves (they are not comparable across families), is not
  available for fits with smooth terms or for split fits, and costs one
  selection plus one joint fit per bootstrap sample.
- `tc_select(split = )` halves the estimation sample; its intervals are
  correspondingly wider. The split is random: another seed can select
  other families.
- Selection of the copula type (Gaussian, Student-t, PCC), of the margins
  or of covariates is not covered by either method.

## Monte Carlo
- `simulation/mc_selection.R` (two scenarios, n = 400, 400 replications,
  bootstrap in 100): the intervals that ignore the selection covered 0.91
  to 0.97 for margin coefficients and Kendall's tau; sample splitting gave
  the same coverage with intervals 42% longer, the bootstrap with selection
  the same coverage within Monte Carlo error. Neither method improved on
  the naive intervals in these scenarios.
# TriCopula 0.13.1

Addendum to run 10: the default covariance for sampling weights, and
Monte Carlo evidence for cases that 0.13.0 left uncovered.

## New
- `tricopula(weights_type = c("auto", "frequency", "sampling"))`. With
  `"auto"` weights that are all whole numbers are frequency weights and
  any others are sampling weights. `fit$weights_type` records the result
  (`NULL` without weights).

## Changed
- **Fits with sampling weights:** `vcov()`, `summary()`, `confint()`,
  `tc_dependence()`, `tc_effects()` and the `model.se` column of
  `tc_bootstrap()` now use the sandwich covariance by default. In 0.13.0
  they used the model-based covariance, whose standard errors were 51 to
  84% of the empirical standard deviation in the Monte Carlo (coverage
  0.69 to 0.91 instead of 0.95). `type = "hessian"` still gives the
  model-based matrix, and `fit$vcov` is unchanged. `summary()` prints the
  kind of weights and the covariance used.
- `type` of these functions now defaults to `NULL` (resolved as above)
  instead of `"hessian"`.

## Unchanged
- Fits without weights, and fits with whole-number weights not declared as
  sampling weights: all results as in 0.13.0 (and, without weights, 0.8.0).
- Estimates and log-likelihoods of all weighted fits.

## Limitations
- Whole-number sampling weights are taken as frequency weights unless
  `weights_type = "sampling"` is given.
- `tc_smooth()` and the tests of smooth terms still use the Bayesian
  covariance, also with sampling weights.
- Missing covariates are still not modelled (rows are dropped); the Monte
  Carlo now shows when that is harmless and when it is not.
# TriCopula 0.13.0

Prior weights and missing responses (extension run 10).

## New
- `tricopula(weights = )`: non-negative prior weights, one per row. The
  log-likelihood is `sum(weights * l_i)`. Integer weights reproduce the
  estimates, log-likelihood and model-based covariance of the data with
  rows repeated. The sandwich covariance uses `sum(w_i^2 s_i s_i')`, i.e.
  it treats each row as one independent unit, as for sampling weights.
  `fit$scores` holds the weighted scores, `fit$loglik_i` the unweighted
  contributions.
- `tricopula(missing = )`: `"fail"` (default, as before), `"complete"`
  (complete-case analysis) or `"likelihood"`: the observed-data likelihood,
  valid when the responses are missing at random. A row with one or two
  missing responses contributes the densities of its observed responses and
  the copula of those responses. Gaussian and Student-t copulas: the
  bivariate margin in closed form. PCC: the pair copula when the root is
  observed, a one-dimensional integral over the root otherwise
  (Gauss-Legendre, 128 nodes mapped to remove the endpoint singularities).
  Rows with a missing covariate or with all responses missing are dropped
  (`fit$rows_dropped`); `fit$observed` is the matrix of observed responses
  and `fit$data` holds the rows used.
- `tc_select(weights = , missing = )`; both bootstraps carry weights and
  missingness patterns.

## Changed
- `.tc_eval()` gains `wt` and `obs`; `.tc_cop_eval()` hands uniforms with
  `NA` to the new `.tc_cop_eval_missing()` (file `R/missing.R`);
  `.tc_fit_margin()` gains weights. None of this is reached with the
  defaults.

## Unchanged
- With `weights = NULL` and `missing = "fail"`: estimates as in 0.8.0 to
  0.12.1 (the 0.8.0 fixtures pass).

## Limitations
- Missing covariates are not modelled: such rows are dropped.
- "Missing at random" is an assumption; nothing in the package tests it.
- A PCC row with a missing root costs about 130 evaluations of the vine
  density, so fits with many such rows are slow. Choosing as root a
  response that is rarely missing avoids the integral.
- Sample averages reported by `tc_dependence()`, `tc_tau()` and
  `fit$correlation` are unweighted.
- `tc_select()` selects families on the complete pairs without weights.
- `tc_diagnostics()`: the Rosenblatt check uses complete rows only.

# TriCopula 0.12.1

Closes limitations of the smooth terms of 0.12.0 (addendum to run 13).

## New
- Smooths with several penalties and with `by` variables: `te()`, `ti()`,
  `s(x, by = g)` (numeric or factor `by`), adaptive and shrinkage bases.
  The generalized Fellner-Schall update uses the pseudo-inverse of each
  smooth's total penalty. `fit$sp` has one entry per penalty,
  `fit$smooth_blocks` one per smooth; `summary()$smooths` has one row per
  smooth with the number of its penalties. `t2()` is still not supported.
- `vcov(fit, type = "corrected")` and `tc_smooth(type = "corrected")`:
  Bayesian covariance plus the first-order correction for
  smoothing-parameter uncertainty of Wood, Pya and Saefken (2016); one
  refit per smoothing parameter.
- `tc_bootstrap(type = "nonparametric")` now works for fits with smooth
  terms: the refits reuse the spline bases of the original fit (new
  argument `smooth_basis` of `tricopula()`). New argument
  `sp = "fixed"` / `"select"` chooses whether the smoothing parameters are
  held fixed or selected again in every refit.
- Smooths without a penalty (`fx = TRUE`) are accepted.

## Changed
- Speed of penalized fits: the inner fits of the smoothing-parameter
  iteration reuse the information matrix of that iteration (Newton steps
  costing one gradient each); the matrix is recomputed only when the
  coefficients have moved; the final covariance uses two Richardson
  iterations instead of four.
- Smoothing-parameter iteration: the log-step of a parameter that keeps
  moving in the same direction is doubled at each iteration (reset when the
  direction changes), and the stopping rule on effective degrees of freedom
  is 0.002 instead of 0.01. Selected smoothing parameters can therefore
  differ slightly from 0.12.0; effects in the null space of the penalty now
  reach it in a few iterations.

## Unchanged
- Formulas without smooth terms: estimates as in 0.8.0 to 0.12.0.

# TriCopula 0.12.0

Smooth terms (extension run 13).

## New
- Penalized smooth terms `s(x)` (mgcv notation) in every formula: means,
  `scale`, `shape`, `dependence` and `df`. Bases and penalties come from
  `mgcv::smoothCon()` with the sum-to-zero constraint absorbed. The
  penalized log-likelihood is maximized with the existing optimizer and the
  smoothing parameters are chosen by the generalized Fellner-Schall
  iteration (Wood and Fasiolo, 2017).
- New argument `sp` of `tricopula()`: fixed smoothing parameters.
- For fits with `s()` terms: `vcov()` is the Bayesian covariance matrix
  `(H + S_lambda)^-1`; `fit$edf` (effective degrees of freedom per
  coefficient), `fit$sp`, `fit$penalties`, `fit$penalized_logLik`,
  `fit$sp_converged`; `logLik()`, `AIC()`, `BIC()` use the total effective
  degrees of freedom; `summary()` prints a table of the smooths with
  effective degrees of freedom and an approximate test.
- `tc_smooth()`: a fitted smooth on a grid with pointwise intervals.
- Unpenalized regression splines (`splines::ns()`, `splines::bs()`) in any
  formula are now covered by tests, including prediction on new data.

## Changed
- `predict()` and the other functions that build design matrices for new
  data now also work when the parametric part of a mean formula has no
  covariates (`y ~ 1`, `y ~ s(x)`).
- `mgcv` added to `Suggests` (needed only for formulas with `s()`).

## Unchanged
- Formulas without `s()` take the unpenalized code path of 0.11.0: the
  0.8.0 regression fixtures pass and the small run of all scripts
  reproduces the 0.8.0, 0.10.0 and 0.11.0 references.

## Limitations
- Only `s()` terms with a single penalty: no `te()`, `ti()`, `t2()`, no
  `by` variables, no random-effect or shrinkage bases with several
  penalties.
- The two-stage (IFM) start and `fit$ifm` are unpenalized.
- Inference for smooths is approximate: the intervals are Bayesian and
  conditional on the smoothing parameters, and the test in `summary()` is a
  Wald statistic referred to a chi-square with the effective degrees of
  freedom.
- `tc_loglik()` and `tc_gradient()` return the unpenalized quantities.
- The nonparametric bootstrap is not available for fits with `s()` terms;
  the parametric bootstrap holds the smoothing parameters fixed.
- `tc_lrtest()` with smooth fits uses the parameter counts, not effective
  degrees of freedom; compare such models by AIC instead.
- Fits with several smooths are slow (tens of seconds at n = 1000).

# TriCopula 0.11.0

Three-parameter margins (extension run 11).

## New
- Margins `"dagum"` (`F(y) = (1 + (y / b)^(-a))^(-p)`) and `"singhmaddala"`
  (`F(y) = 1 - (1 + (y / b)^a)^(-q)`): scale `b = exp(eta)`, `a = exp(.)`
  through `scale`, and a third parameter `p` or `q` on the log scale with
  its own linear predictor through the new argument `shape = ~ ...`. All
  scores are analytic. They coincide with GJRM's `DAGUM` and `SM`.
- Margin `"gpareto"`: generalized Pareto for exceedances `y > 0`, scale
  `sigma = exp(eta)`, shape `xi = exp(.) - 0.5 > -0.5` through `scale`
  (two parameters; analytic scores). It coincides with GJRM's `GPII`.
- `tricopula(shape = )` and `tc_select(shape = )`: one-sided formula or a
  list of three for the third predictor. Coefficient names
  `margin1_log_p`, `margin1_log_q` (with covariates `margin1_log_p_x`).
  `tc_dependence()` and the bootstrap report the third parameters
  (`margin1_p`, `margin2_q`) after the scale/shape rows; `predict()`,
  `simulate()`, `tc_jointprob()` and `tc_effects()` use them.

## Changed
- The margins have two or three predictors. `.tc_index()` gains `shapel`
  (positions of the third predictor inside each margin block, empty for
  two-parameter margins); `.tc_margin_eval()`, `.tc_margin_cdf()`,
  `.tc_margin_quantile()`, `.tc_margin_mean()` and `.tc_margin_pars()` gain
  an optional `shape` argument. The parameter layout of models with
  two-parameter margins is unchanged.
- `VGAM` added to `Suggests` (reference densities in the tests).

## Unchanged
- With two-parameter margins the code path, names and estimates are those
  of 0.10.0: the 0.8.0 regression fixtures pass and the small run of all
  scripts reproduces the 0.8.0 and 0.10.0 references.

## Limitations
- No threshold argument for the generalized Pareto margin: the response
  must be the exceedance over a threshold chosen by the user. The shape is
  restricted to `xi > -0.5`.
- The means of the new margins are infinite for part of the parameter space
  (`a <= 1` Dagum, `a q <= 1` Singh-Maddala, `xi >= 1` generalized Pareto);
  `predict(type = "mean")` then returns `Inf`.
- The two shape parameters of the Dagum and Singh-Maddala margins are
  weakly identified in small samples (see the Monte Carlo of run 11).

# TriCopula 0.10.0

Copula extensions (extension run 12).

## New
- Pair-copula families `"frank"` (parameter `theta` on the identity scale,
  `|theta| <= 100`, negative values give negative dependence) and `"joe"`
  with rotations `"joe90"`, `"joe180"`, `"joe270"` (`log(theta - 1)`):
  log-density, h-function, Kendall's tau, tail dependence, sampling. The
  Frank Kendall's tau uses the Debye function by Gauss-Legendre quadrature
  (agreement with `stats::integrate` below 1e-10).
- `tricopula(df = ~ ...)`: linear predictor for `log(df - 2.05)` of the
  Student-t copula, or of each Student-t pair of a PCC (a list of three
  formulas gives one per pair). `fit$df` is then the sample average and
  `fit$df_i` holds the values per observation; `predict(type = "dependence")`
  gains `df` columns; `simulate()`, `tc_jointprob()`, `tc_effects()`,
  `tc_dependence()`, `tc_diagnostics()` and `tc_bootstrap()` use the
  observation-specific values.
- `tricopula(root = )`: root of the C-vine (1, 2 or 3). The pairs are
  (root, a), (root, b), (a, b | root) with a < b; coefficient names, pair
  labels and all outputs follow the vine (e.g. `pcc13_given_2_...`,
  `tau13|2`). `tc_rpcc(root = )` simulates from such a vine.
- `tc_select(root = )`: selection for a given root, or `root = "auto"`,
  which selects families and fits jointly for each of the three roots and
  returns the root with the smallest AIC/BIC (three joint fits).
  `candidates = "all"` adds Frank and Joe to the default candidates. The
  returned list gains the element `root` (and `roots` for `"auto"`).

## Changed
- `.tc_index()` gains the list `dfl` (positions of each df predictor). The
  parameter layout is unchanged; df coefficients stay at the end, and
  intercept-only df predictors keep their names (`log_df_minus_2.05`,
  `pcc12_log_df_minus_2.05`).
- `"frank"` is now a valid family, so the three tests that used it as the
  example of an invalid family name use `"amh"` instead.

## Unchanged
- With the defaults (`df = ~ 1`, `root = 1`, default `candidates`) the code
  path, the parameter names and the estimates are those of 0.9.0: the 0.8.0
  regression fixtures pass, and the small run of all scripts reproduces the
  0.8.0 reference.

## Limitations
- Derivatives of the Frank and Joe pair copulas use the same central
  differences as the other pair copulas.
- `tc_taildep()` and `tc_tau()` use sample-average predictors, also for a
  covariate-dependent df.
- The Frank inverse h-function loses accuracy where the h-function is flat
  in double precision (very strong dependence, `|theta|` above about 30).
- Standard errors after `tc_select()` ignore the selection of families and
  root.

# TriCopula 0.9.0

Baseline and housekeeping release. No file in `R/` was changed.

## New
- Numeric regression fixtures: `tests/testthat/fixtures/regression_0.8.0.rds`
  (six data sets covering the Gaussian and Student-t copulas, the
  Weibull / inverse Gaussian / Beta margins, `scale` and `dependence`
  formulas, and two pair-copula constructions) with the results of 0.8.0,
  and `tests/testthat/test-regression.R`, which requires the default fit to
  reproduce them to 1e-6 and the likelihood at a fixed parameter vector to
  1e-10.
- GitHub Actions workflow `.github/workflows/R-CMD-check.yaml` (Windows,
  macOS, Linux; R release, devel, oldrel). Not yet executed.

## Changed
- `DESCRIPTION`: author and maintainer Ibrahim Niankara. `LICENSE`:
  copyright holder Ibrahim Niankara (licence still MIT).
- `.gitignore` and `.Rbuildignore` extended for a git repository rooted at
  the package directory.

## Unchanged
- Model, likelihood and estimates are identical to 0.8.0.

# TriCopula 0.8.0

## New
- `tc_select()`: sequential selection of the three pair-copula families of a
  PCC by AIC or BIC, followed by a joint fit of the selected model.
- Vignette "Trivariate copula regression with TriCopula": a worked example
  on the 1990 cross-section of the Spanish manufacturing firms in
  `plm::Snmesp` (output, employment and wages given capital).
- The package now passes the full `R CMD check`, including the PDF manual
  and the vignette build (the vignette needs pandoc).
- New suggests: `knitr`, `rmarkdown`, `plm` (vignette only).

## Unchanged
- Model, likelihood and estimates are identical to 0.7.0.

# TriCopula 0.7.0

## New
- Pair-copula construction: `tricopula(..., copula = "pcc", families = c(f12,
  f13, f23|1))`, a simplified C-vine rooted at the first response, estimated
  jointly with the margins. Families: `gaussian`, `student`, `clayton`,
  `gumbel` and the 90/180/270-degree rotations of Clayton and Gumbel.
  `dependence` covariates, `scale` covariates, bootstrap, effects, joint
  probabilities, simulation and diagnostics all work with a PCC.
- `tc_rpcc()`: simulate from a trivariate C-vine.
- For a PCC, `tc_tau()` returns the three pair taus, `tc_taildep()` the
  lower/upper coefficients of the two unconditional pairs, and
  `tc_diagnostics()` the largest correlation among the Rosenblatt
  transforms.
- New suggest: `VineCopula` (tests only).

## Fixed
- Spurious very large likelihoods at near-singular correlation matrices.
  Since 0.3.0 the copula density used a closed-form 3 x 3 inverse whose
  determinant cancels catastrophically when a correlation is numerically
  +/-1; the density then evaluated to rounding noise and the optimizer
  could be attracted to it (2 of 300 Student-t fits in the 0.6.0 study, and
  the unreliable `fit$ifm$logLik` seen when the df ran to infinity). The
  determinant now uses the product form
  `(1 - rho12^2)(1 - rho13^2)(1 - rho23|1^2)`, and parameter values with
  `|atanh(rho)| > 6` (`|rho| > 0.999988`) are rejected. The same bound
  applies to the predictors of PCC pairs.

## Unchanged
- Gaussian and Student-t copula fits are identical to 0.6.0. An all-Gaussian
  PCC reproduces the Gaussian-copula fit.

# TriCopula 0.6.0

## New
- Three more margin families, each usable in any equation and with `scale`
  covariates: `"weibull"` (scale `exp(eta)`, shape), `"invgauss"` (mean
  `exp(eta)`, sigma; variance `mu^3 sigma^2`) and `"beta"` for responses in
  (0, 1) (mean `plogis(eta)`, precision). Parameterizations agree with
  GJRM's `WEI` and `IG` families; GJRM's `BE` uses `sigma` with
  `precision = (1 - sigma^2) / sigma^2`.
- `tc_effects()`: average marginal effects of a covariate on the three
  conditional means or on a joint probability, as an average derivative or
  a discrete change, with delta-method standard errors.
- `predict(type = "mean")` and `predict()` cover the new families.

## Internal
- Family-specific code is gathered in `R/margins.R` (density, CDF, quantile,
  mean, scores, starting values, labels). For the inverse Gaussian and Beta
  margins the CDF derivatives use central differences.

## Unchanged
- Existing families, parameter names and estimates are identical to 0.5.0.

# TriCopula 0.5.0

## New
- `tc_bootstrap()`: parametric (simulate from the fit) and nonparametric
  (resample rows) bootstrap standard errors and percentile intervals for all
  coefficients and natural-scale parameters, optionally on a `parallel`
  cluster. Refits start from the original estimates.
- The fit object now stores `data` and `control` (needed for refitting).
- New import: `parallel` (base R).

## Fixed
- Student-t copula fitted to data without tail dependence: the df drifts
  towards infinity, where the likelihood is flat, and 0.4.0 then returned no
  covariance matrix. The fit now warns that it is practically Gaussian,
  returns the covariance of all other parameters (conditional on df) and
  reports `NA` for the df standard error; `fit$df_boundary` flags the case.

## Unchanged
- Model, likelihood, parameterization and all estimates are identical to
  0.4.0.

# TriCopula 0.4.0

## New
- Covariates on the marginal scale and shape: `tricopula(..., scale = ~ x)`
  (or a list of three one-sided formulas) gives `log sigma` (Gaussian,
  log-normal) or `log shape` (Gamma) of each margin its own linear
  predictor, i.e. a location-scale regression in every margin. Analytic
  scores cover the new coefficients.
- `predict()`, `simulate()` and `tc_jointprob()` use row-specific scale and
  shape.

## Unchanged
- The default `scale = ~ 1` reproduces 0.3.0 exactly, including parameter
  names and order. With covariates the coefficients are named
  `margin1_log_scale_<term>` / `margin3_log_shape_<term>`.
- With covariates on the scale, `tc_dependence()` and `summary()` report the
  sample-average scale or shape.

# TriCopula 0.3.0

## New
- Covariate-dependent dependence: `tricopula(..., dependence = ~ x)` (or a
  list of three one-sided formulas) gives each of `atanh(rho12)`,
  `atanh(rho13)` and `atanh(rho23|1)` its own linear predictor. The
  correlation matrix is observation-specific and positive definite for every
  observation. Analytic scores cover the new coefficients.
- `predict(type = "dependence")`: fitted correlations, partial correlation
  and Kendall's tau per covariate row. `fit$rho` holds them for the sample.
- `tc_lrtest()`: likelihood-ratio test of nested fits, e.g. constant against
  covariate-dependent dependence.
- `simulate()` and `tc_jointprob()` honour row-specific dependence.

## Unchanged
- The default `dependence = ~ 1` reproduces 0.2.0 exactly, including
  parameter names and order. With covariates the dependence coefficients are
  named `atanh_partial_rho12_<term>` etc.
- With covariate-dependent dependence, `fit$correlation`, `tc_dependence()`,
  `tc_tau()` and `tc_taildep()` report sample-average correlations.

## Internal
- The copula density uses the closed-form 3 x 3 inverse and determinant so
  that it is vectorized over observation-specific correlations.
- The two-stage (IFM) copula step now uses the analytic copula score.
- The Student-t degrees of freedom remain constant.

# TriCopula 0.2.0

First executed and tested version (R 4.6.2, Windows). The model, margins,
copulas and parameterization of 0.1.0 are unchanged.

## New
- Analytic score for the full joint log-likelihood (central differences only
  for the Gamma CDF shape derivative and the Student-t degrees of freedom);
  verified against numerical differentiation.
- Two-stage (IFM) starting values: marginal ML per equation, then copula ML.
  The IFM estimates and log-likelihood are returned in `fit$ifm`.
- `nstart` multi-start optimization and Newton refinement (`polish`).
- Observed information from the Jacobian of the analytic gradient;
  `vcov(type = "sandwich")` robust covariance; `confint()`.
- `tc_dependence()`: correlations, partial correlation, degrees of freedom
  and marginal scale/shape with delta-method standard errors and intervals.
- `tc_diagnostics()`, `tc_gradient()`, `tc_tau()`, `tc_taildep()`,
  `tc_jointprob()`, `tc_compare()`, `tc_rcopula()`, `simulate()`,
  `predict(type = "mean")`, `nobs()`.
- Help pages for every exported function; `R CMD check` passes.

## Changed relative to 0.1.0
- `tricopula()` gains `nstart` and `polish`; `predict()` gains a default
  `newdata = NULL` and `type`; `vcov()` gains `type`.
- Default starting values are now the IFM estimates instead of OLS
  coefficients with zero dependence parameters.
- The returned `hessian` is computed with `numDeriv` from the analytic
  gradient rather than by `optim()`.
- New imports: `numDeriv`. New suggests: `mvtnorm`.

## Fixed
- `predict.tricopula()` in 0.1.0 referenced an object defined inside a
  nested block and ignored factor levels of the estimation sample.

# TriCopula 0.1.0

- Initial research prototype created (never executed).
- Three continuous regression margins: Gaussian, log-normal, Gamma.
- Trivariate Gaussian and Student-t copula likelihoods.
- Positive-definite correlation parameterization using partial correlations.
- BFGS estimation, numerical Hessian, basic prediction and PIT extraction.
- Simulation helper and initial Monte Carlo template.
