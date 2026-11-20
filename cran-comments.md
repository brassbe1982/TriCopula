## Submission

This is a first submission of TriCopula (version 0.17.3).

## Test environments

* Local: Windows 11 x64 (build 26200), R 4.6.2 (2026-10-29 ucrt),
  `R CMD check --as-cran` on the tarball built by `R CMD build`.

No other platform has been checked yet. [Maintainer to add before
submission: the results of win-builder (R-devel), of R-hub or the GitHub
Actions workflow in `.github/workflows/` (Linux, macOS), and of the macOS
builder, once these have been run.]

## R CMD check results

0 errors | 0 warnings | 2 notes

* `checking CRAN incoming feasibility ... NOTE`
  `Maintainer: 'Ibrahim Niankara <Ibrahim.niankara@aaua.ac.ae>'`
  `New submission`

  This is a first submission.

* `checking HTML version of manual ... NOTE`
  `Skipping checking math rendering: package 'V8' unavailable`

  V8 is not installed on the local machine; the note is about the check
  environment, not the package.

Words in DESCRIPTION that a spell checker may flag are proper names
(Aas, Czado, Frigessi, Bakken, Marra, Radice, Clayton, Gumbel, Frank, Joe,
Pareto, Dagum, Singh, Maddala) or terms of the field (copula, C-vine).

## Timing

The complete check took 2.7 minutes locally: tests 58 seconds, vignette
30 seconds, examples a few seconds. The 47 slowest tests (about 16 minutes
in total) call `testthat::skip_on_cran()`; they run when `NOT_CRAN` is
`"true"` and pass (107 tests, 1779 expectations).

## Other

* No compiled code. Imports only base and recommended packages and
  numDeriv; all packages in Suggests are used conditionally
  (`requireNamespace()` in the code, `skip_if_not_installed()` in the tests,
  and the vignette is evaluated only if plm is available). The tests were
  run in CRAN mode against a library holding only the hard dependencies
  and testthat (with mgcv and survival also hidden in a second run): no
  failures or errors.
* Examples, tests and the vignette use one core; parallel computation
  happens only on a cluster supplied by the user (`tc_bootstrap(cl = )`).
* Nothing is written outside `tempdir()`, and the vignette restores the
  options it sets.
* [Maintainer to add before submission: `URL` and `BugReports` fields in
  DESCRIPTION, once the public repository exists.]
