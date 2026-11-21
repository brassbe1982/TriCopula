# TriCopula 0.17.5 - research prototype

Full maximum-likelihood regression for three continuous responses linked by
a trivariate Gaussian or Student-t copula, or by a pair-copula construction
with Gaussian, Student-t, Clayton, Gumbel, Frank and Joe pair copulas.

## Model

For observation `i` and responses `j = 1, 2, 3`:

    log f(y1, y2, y3 | x) = sum_j log f_j(y_j | x) + log c(u1, u2, u3),   u_j = F_j(y_j | x)

| Margin      | Model                                   | Extra parameter (log scale) |
|-------------|-----------------------------------------|-----------------------------|
| `gaussian`  | `y ~ N(x'b, sigma^2)`                   | `log sigma`                 |
| `lognormal` | `log y ~ N(x'b, sigma^2)`               | `log sigma`                 |
| `gamma`     | `E(y) = exp(x'b)`, shape `a`            | `log a`                     |
| `weibull`   | scale `exp(x'b)`, shape `k`             | `log k`                     |
| `invgauss`  | `E(y) = exp(x'b)`, `Var = mu^3 sigma^2`   | `log sigma`                 |
| `beta`      | `E(y) = plogis(x'b)`, precision `phi`   | `log phi`                   |
| `gpareto`   | scale `exp(x'b)`, shape `xi > -0.5`, `y > 0` an exceedance | `log(xi + 0.5)` |
| `dagum`     | `F(y) = (1 + (y/b)^(-a))^(-p)`, `b = exp(x'b)` | `log a`; third: `log p` |
| `singhmaddala` | `F(y) = 1 - (1 + (y/b)^a)^(-q)`, `b = exp(x'b)` | `log a`; third: `log q` |

The Dagum and Singh-Maddala margins have a third parameter with its own
linear predictor, `shape = ~ x` (or a list of three formulas; default `~ 1`).
These three margins use the parameterizations of GJRM's `GPII`, `DAGUM` and
`SM`.

The three margins may be any mixture of these. The extra parameter of each
margin can have its own linear predictor, `log sigma_i = z_i' g` (or
`log a_i` for the Gamma margin), specified with `scale = ~ x` or a list of
three one-sided formulas; the default `~ 1` keeps it constant. The copula correlation matrix
is parameterized by `atanh(rho12)`, `atanh(rho13)` and `atanh(rho23|1)`, with

    rho23 = rho12 rho13 + sqrt((1 - rho12^2)(1 - rho13^2)) rho23|1,

so every real parameter vector gives a positive-definite matrix
(`|R| = (1 - rho12^2)(1 - rho13^2)(1 - rho23|1^2) > 0`). The Student-t
degrees of freedom are `2.05 + exp(.)`; `df = ~ x` gives them a linear
predictor, `log(df_i - 2.05) = v_i' z` (default `~ 1`, constant).

Each dependence parameter can have its own linear predictor,

    atanh(rho12_i) = w_i' d12,  atanh(rho13_i) = w_i' d13,  atanh(rho23|1_i) = w_i' d23,

specified with `dependence = ~ x` (one formula for all three) or a list of
three one-sided formulas. The default `~ 1` is constant dependence.

## Smooth terms

Every formula (means, `scale`, `shape`, `dependence`, `df`) may contain
penalized smooths in mgcv notation, for example
`y1 ~ z + s(x)` or `dependence = ~ s(x, k = 6)`. The bases and penalties
come from `mgcv::smoothCon()`; the penalized log-likelihood
`l - 1/2 sum_k lambda_k b_k' S_k b_k` is maximized and the smoothing
parameters are chosen by the generalized Fellner-Schall iteration (Wood and
Fasiolo, 2017), or fixed with `sp`. For such fits `vcov()` is the Bayesian
covariance matrix, `AIC()` uses the effective degrees of freedom,
`summary()` lists the smooths, and `tc_smooth()` returns a fitted curve with
pointwise intervals. Supported: `s()`, `te()`, `ti()`, with or without a
`by` variable, including smooths with several penalties (not `t2()`).
`vcov(fit, type = "corrected")` adds a first-order correction for the
uncertainty of the smoothing parameters, and both bootstraps work for such
fits. Unpenalized splines (`splines::ns(x, 4)`) work in every formula
without mgcv.

## Weights and missing responses

`weights` gives non-negative prior weights (frequency or sampling weights):
the log-likelihood is the weighted sum of the rows' contributions. Weights
that are not whole numbers are taken as sampling weights, and `summary()`,
`vcov()`, `confint()`, `tc_dependence()` and `tc_effects()` then use the
sandwich covariance by default (0.13.1); whole-number weights are taken as
frequency weights with the model-based covariance. `weights_type` overrides
this.

`missing = "likelihood"` keeps rows in which one or two of the three
responses are missing and uses the observed-data likelihood: the densities
of the observed responses and the copula of those responses. This is valid
when responses are missing at random, where complete-case analysis
(`missing = "complete"`) can be biased. Rows with missing covariates are
dropped. In a PCC a missing root is integrated out numerically, which is
slow; choose as root a response that is rarely missing.

## Pair-copula construction

`copula = "pcc"` with `families = c(f12, f13, f23|1)` uses the simplified
C-vine rooted at the first response,

    c(u1, u2, u3) = c12(u1, u2) c13(u1, u3) c23|1(h12(u2 | u1), h13(u3 | u1)),

with each pair one of `gaussian`, `student`, `clayton`, `gumbel` or a
90/180/270-degree rotation of Clayton or Gumbel. Clayton gives lower-tail and
Gumbel upper-tail dependence; 180-degree rotations swap the tails and 90/270
give negative dependence. The pair parameters use the links `atanh(rho)`,
`log(theta)` (Clayton) and `log(theta - 1)` (Gumbel), and can depend on
covariates through `dependence`. Three Gaussian pairs reproduce the
trivariate Gaussian copula.

Since 0.10.0 the pairs may also be `frank` (parameter `theta` itself; any
sign; no tail dependence) and `joe` with rotations `joe90`, `joe180`,
`joe270` (`log(theta - 1)`; upper-tail dependence). `df = ~ x` (or a list of
three formulas) gives each Student-t pair a predictor for
`log(df - 2.05)`.

The root of the vine is the first response unless `root = 2` or `root = 3`
is given. With root `r` the pairs are (r, a), (r, b) and (a, b | r), with
`a < b` the other two responses, in this order in `families`, `dependence`
and the coefficient names (`pcc12_...`, `pcc23_...`, `pcc13_given_2_...` for
`root = 2`). The root is the first argument of the two unconditional pair
copulas, which only matters for 90/270-degree rotations.
`tc_select(root = "auto")` selects families for each root, fits the three
models jointly and returns the one with the smallest AIC or BIC.

Since 0.15.0 the vine need not be simplified. With `simplified = FALSE` the
predictor of the conditional pair is `w' delta + f(u_r)`, with `f` a
penalized cubic B-spline in the probability integral transform of the root
(`k_cond` basis functions, 10 by default since 0.17.0, fixed knots, `f`
integrating to zero). The
margins, the three pairs and `f` are estimated jointly; `u_r` moves with the
parameters of the root's margin. `tc_smooth()` returns the fitted `f`,
`tc_simplifying_test()` tests the simplified vine against a spline
alternative by a likelihood-ratio test.

Since 0.16.0 responses may be truncated and censored. `truncation` gives
known bounds per response (numbers or columns of the data). With
`truncation_type = "margin"` the margin is the truncated distribution and
the model describes the truncated population; with `"joint"` the
untruncated model holds and rows are observed only if the truncated
response is inside its bounds (one truncated response). `censoring` names a
status column per response, coded as in `survival::Surv()` (1 exact, 0
right-, 2 left-, 3 interval-censored). A row contributes the densities of
its exact responses and the copula density integrated over its censored
ones; the integral is computed by Gauss-Legendre quadrature (`cens_nodes`)
for every copula and every pattern of censored responses. Since 0.17.0
one censored response per row is integrated in closed form
(`cens_method = "closed"`, the default since 0.17.1; `"quadrature"`
reproduces 0.16.0), and every fit with censoring checks its quadrature against
twice the nodes (`tc_cens_check()`, `fit$cens_accuracy`) and warns if an
estimate would move by more than 0.1 standard errors.

## Estimation

1. Two-stage (IFM) start: marginal ML per equation, then copula ML with the
   margins fixed. Kept in `fit$ifm` for comparison.
2. Joint BFGS maximization of the full log-likelihood over all marginal and
   dependence parameters with an analytic gradient, optionally from several
   starts (`nstart`), followed by Newton refinement.
3. Observed information = Jacobian of the analytic gradient; covariance =
   its inverse, or the sandwich `H^-1 S'S H^-1`.

## Installation and use

Source and issue tracker: <https://github.com/brassbe1982/TriCopula>.

```r
# from GitHub
install.packages("remotes")
remotes::install_github("brassbe1982/TriCopula")
# or from a source tarball:  R CMD INSTALL TriCopula_0.17.5.tar.gz
install.packages("mvtnorm")   # suggested; needed for tc_jointprob()

library(TriCopula)
fit <- tricopula(
  formulas = list(income ~ education + age,
                  consumption ~ education + age,
                  saving ~ education + age),
  data = d,
  margins = c("lognormal", "gamma", "gaussian"),
  copula = "student", nstart = 2,
  scale = ~ education,           # optional; default ~ 1
  dependence = ~ education)      # optional; default ~ 1

summary(fit)                 # coefficients + natural-scale dependence table
tc_diagnostics(fit)          # convergence, gradient, information, PIT checks
tc_dependence(fit)           # rho, rho23|1, df with delta-method intervals
vcov(fit, type = "sandwich") # robust covariance
tc_tau(fit); tc_taildep(fit)
predict(fit, newdata, type = "mean")
tc_jointprob(fit, thresholds = c(4, 2.5, 0), newdata)
simulate(fit, seed = 1)
tc_compare(fit_gaussian, fit)
predict(fit, newdata, type = "dependence")  # rho and tau per row
tc_lrtest(fit_constant, fit)                # H0: constant dependence
tc_bootstrap(fit, B = 200, type = "nonparametric")  # bootstrap SEs and intervals
tc_effects(fit, "education")                # average marginal effects on the means

# asymmetric tail dependence: pair-copula construction
fit_pcc <- tricopula(formulas, d, margins = c("lognormal", "gamma", "gaussian"),
                     copula = "pcc",
                     families = c("clayton", "gumbel180", "gaussian"))
tc_compare(fit, fit_pcc); tc_taildep(fit_pcc)

# or let the data choose the three pair families (sequentially, by AIC)
sel <- tc_select(formulas, d, margins = c("lognormal", "gamma", "gaussian"))
sel$families; summary(sel$fit)

# 0.10.0: Frank and Joe pairs, another root, covariate-dependent df
tricopula(formulas, d, margins = c("lognormal", "gamma", "gaussian"),
          copula = "pcc", families = c("frank", "joe180", "student"),
          root = 2, df = ~ education)
tricopula(formulas, d, margins = c("lognormal", "gamma", "gaussian"),
          copula = "student", df = ~ education)
tc_select(formulas, d, margins = c("lognormal", "gamma", "gaussian"),
          candidates = "all", root = "auto")$roots

# 0.11.0: three-parameter margins; the Dagum p depends on education
tricopula(formulas, d, margins = c("dagum", "singhmaddala", "gaussian"),
          shape = list(~ education, ~ 1, ~ 1))

# 0.12.0: penalized smooths in the mean and in the dependence
fs <- tricopula(list(income ~ education + s(age), consumption ~ education + age,
                     saving ~ education + age), d,
                margins = c("lognormal", "gamma", "gaussian"),
                dependence = list(~ s(age), ~ 1, ~ 1))
summary(fs)$smooths; AIC(fs)
plot(fit ~ age, tc_smooth(fs, "margin1:s(age)"), type = "l")

# 0.13.0: survey weights, and responses that are partly missing
fw <- tricopula(formulas, d, margins = c("lognormal", "gamma", "gaussian"),
                weights = d$wgt, missing = "likelihood")
fw$rows_dropped; colMeans(fw$observed); fw$weights_type; summary(fw)

# 0.14.0: inference that allows for the selection of the families
sel <- tc_select(formulas, d, margins = c("lognormal", "gamma", "gaussian"))
tc_bootstrap(sel$fit, B = 200, type = "nonparametric", select = TRUE)
sp <- tc_select(formulas, d, margins = c("lognormal", "gamma", "gaussian"),
                split = 0.5, seed = 1)   # select on one half, estimate on the other
summary(sp$fit)

# 0.14.1: margins, covariates and copula type selected as well; tail
# dependence, means and joint probabilities after selection
sa <- tc_select(formulas, d, margins = list(c("lognormal", "gamma"), "gamma", "gaussian"),
                covariates = "backward", copula = "auto")
ba <- tc_bootstrap(sa$fit, B = 200, type = "nonparametric", select = TRUE,
                   thresholds = c(1, 1, 0), newdata = d[1:5, ])
ba$derived; ba$copulas; ba$margin_selection; ba$covariate_selection

# 0.15.0: non-simplified vine and the test of the simplifying assumption
ns <- tricopula(formulas, d, margins = c("lognormal", "gamma", "gaussian"),
                copula = "pcc", families = sel$families, simplified = FALSE)
summary(ns)$smooths; tc_smooth(ns, 1); tc_simplifying_test(sel$fit)

# 0.16.0: censoring (status 1 exact, 0 right-censored) and truncation
# d$s1 is the status of y1; rows are in the sample only if y3 > 0.5
fc <- tricopula(formulas, d, margins = c("lognormal", "gamma", "gaussian"),
                censoring = list(c(status = "s1"), NULL, NULL),
                truncation = list(NULL, NULL, c(0.5, Inf)), truncation_type = "joint")

# 0.17.0: accuracy of the quadrature over censored responses
fc$cens_accuracy; tc_cens_check(fc, nodes = c(200, 100, 40))
```

`vignette("TriCopula")` works through an empirical example (Spanish
manufacturing firms: output, employment and wages given capital).

## Exported functions

`tricopula`, `tc_loglik`, `tc_gradient`, `tc_dependence`, `tc_diagnostics`,
`tc_tau`, `tc_taildep`, `tc_jointprob`, `tc_compare`, `tc_pit`,
`tc_joint_logscore`, `tc_lrtest`, `tc_bootstrap`, `tc_effects`, `tc_select`, `tc_smooth`,
`tc_simplifying_test`, `tc_rcopula`, `tc_rpcc`,
`tc_simulate`; methods `print`,
`summary`, `logLik`, `nobs`, `vcov`, `confint`, `predict`, `simulate`.

## Relation to GJRM

Established (GJRM 0.2-6.9, checked against the installed package and its
source): `gjrm.pcc()` takes three fitted `GJRM::gamlss()` objects, holds
their PITs fixed, and estimates three bivariate pair copulas one after the
other (`gjrm(..., fixed.margins = TRUE)`), with family selection, optional
covariate effects in the copula parameter and, unless `simplified = TRUE`, a
smooth of the conditioning PIT in the third copula. It is a sequential
estimator with far richer margins and copula families than this package.

TriCopula's distinct feature is narrower: all marginal and dependence
parameters are estimated jointly from one likelihood, with a joint
covariance matrix, and the Student-t degrees of freedom are estimated rather
than supplied or chosen from a grid.

For the Frank and Joe families and for a vine rooted at another variable
(GJRM: `copulas = c("F", "J0", ...)`, `cond.var`), the two-stage estimates of
GJRM evaluated in TriCopula's likelihood reproduce the log-likelihood GJRM
reports to 1e-12 (`validation/benchmark_gjrm_pcc_frankjoe.R`). A linear
predictor for the degrees of freedom has no counterpart in `gjrm.pcc()`,
which takes them as given or from a grid.

When the vine is simplified, rooted at variable 1 and uses Gaussian pair
copulas, both packages target the same model, with or without a linear
predictor for the copula parameters (`eq.theta` in GJRM, `dependence` here).
`validation/benchmark_gjrm_pcc.R`, `benchmark_gjrm_pcc_theta.R` and
`benchmark_gjrm_pcc_scale.R` (the last with sigma equations in the margins)
compare them parameter by parameter.

The non-simplified vine (0.15.0) is established GJRM functionality: it is
the default of `gjrm.pcc()`, which adds `s(u)` of the conditioning PIT to the
parameter equation of the third copula, with the PIT taken from the fitted
margin. TriCopula adds its joint estimation with the margins (the PIT then
depends on the parameters being estimated) and a likelihood-ratio test of
the simplifying assumption. `validation/benchmark_gjrm_nonsimplified.R`
compares the two; the copula log-likelihoods GJRM reports are reproduced
with TriCopula's pair copulas and GJRM's per-observation parameters.

Censoring (0.16.0): GJRM handles right, left, interval and mixed censoring
for its survival margins (link-based models with a monotone spline of
time, `family = "-cloglog"`, `"-logit"`, `"-probit"`), and left truncation
through `truncation.time`; `left.trunc` is for count distributions. These
are not the parametric margins of TriCopula, so the two cannot be compared
parameter by parameter. TriCopula's censored likelihood is checked instead
against `survival::survreg()` (marginal stage) and against direct numerical
integration with `mvtnorm` and `VineCopula`.

## Limitations

1. Smooth terms: no `t2()`. Inference for smooths is approximate: Bayesian
   intervals, optionally with a first-order correction for
   smoothing-parameter uncertainty, and an approximate test. Fits with
   several smooths are several times slower than GJRM's sequential fit.
2. With covariate-dependent dependence or df the summary functions report
   sample averages. A df predictor needs large samples: the df is weakly
   identified, most of all where it is large.
3. Trivariate Gaussian and Student-t copulas, and a C-vine with Gaussian,
   Student-t, Clayton, Gumbel, Frank and Joe pairs, simplified or not. In
   the non-simplified vine only the dependence parameter of the conditional
   pair depends on the conditioning value. `tc_select()` chooses
   families sequentially and, on request, the root, the margins, the terms
   of the mean formulas and the copula type by AIC/BIC; the standard
   errors of that fit do not account for these choices. `tc_select(split = )`
   and `tc_bootstrap(select = TRUE)` do (0.14.0, 0.14.1); the bootstrap gives
   no interval after selection for the copula parameters themselves.
4. Nine margin families: Gaussian, log-normal, Gamma, Weibull, inverse
   Gaussian, Beta, generalized Pareto (two parameters), Dagum and
   Singh-Maddala (three parameters). The generalized Pareto margin has no
   threshold argument (supply exceedances) and its shape is restricted to
   `xi > -0.5`. The two shape parameters of the Dagum and Singh-Maddala
   margins are weakly identified in small samples.
5. Censoring and truncation (0.16.0) use numerical integration over the
   censored responses: rows with two or three censored responses are slow,
   and the accuracy drops under very strong dependence (see `cens_nodes`;
   since 0.17.0 the fit checks this and warns, and one dimension of the
   integration is done in closed form). In a vine there is no closed
   form for a censored root.
   Censoring is assumed independent of the responses given the covariates.
   Truncation of type `"joint"` handles one truncated response. Missing
   covariates are not modelled (rows dropped); the observed-data likelihood
   for missing responses assumes missing at random and cannot be combined
   with censoring.
6. BFGS works on the raw covariate scale; badly scaled covariates should be
   standardized by the user.
7. PITs are clipped to `[1e-12, 1 - 1e-12]` before inverse-CDF evaluation.
8. `tc_jointprob()` evaluates one trivariate integral per distinct covariate
   row and is slow for large `newdata` with a Student-t copula.
9. The KS tests in `tc_diagnostics()` ignore parameter estimation.
10. `R CMD check` passes on Linux, macOS and Windows with R release and
    oldrel, and on Linux and Windows with R-devel (GitHub Actions,
    `.github/workflows/R-CMD-check.yaml`). Standard errors come from a
    numerically differentiated score and agree across platforms to about
    3e-4 relative; estimates agree to 1e-6 or better.

## Acknowledgment

This R package is brought to the world freely and open-access through the
generous funding of BRASS Digital Lab, an LLC registered in Burkina Faso
under Trade Register (RCCM) BF-OUA-01-2026-B13-02466 and Social Security
Fund (CNSS) 1493318 R, with registered office at 02 BP 5444, Ouagadougou 02,
Secteur 33, Lot 13, Section 762, Lot 28, Burkina Faso. Website:
<https://brassbe1982-bdlweb.share.connect.posit.cloud/>.

## Regression baseline

`tests/testthat/fixtures/regression_0.8.0.rds` holds six small data sets
with the coefficients, log-likelihoods, standard errors and fixed-point
likelihood and score values produced by version 0.8.0.
`tests/testthat/test-regression.R` requires later versions to reproduce them
(fits to 1e-6, fixed-point likelihood to 1e-10). Version 0.9.0 changes no
model code: its estimates are those of 0.8.0. Versions 0.10.0 and 0.11.0 add
arguments whose defaults (`df = ~ 1`, `root = 1`, `shape = ~ 1`) give the
same estimates.
