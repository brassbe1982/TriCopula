## Submission

This is a first submission of TriCopula (version 0.17.4).

## Test environments

* Local: Windows 11 x64 (build 26200), R 4.6.2 (2026-10-29 ucrt),
  `R CMD check --as-cran` on the tarball built by `R CMD build`.
* GitHub Actions (`R CMD check --as-cran --no-manual`, complete test suite
  with `NOT_CRAN = "true"`), run 2 of
  <https://github.com/brassbe1982/TriCopula/actions>:
  - Ubuntu 26.04, R-devel (2026-11-19 r90939), R 4.6.2 and R 4.5.3: OK
  - macOS (arm64), R 4.6.2 and R 4.5.3: OK
  - Windows Server 2022, R-devel, R 4.6.2 and R 4.5.3: OK
  - macOS, R-devel: not checked; the runner could not compile the
    dependency Matrix (`libintl.h` not found) before the check started.
* win-builder, R-devel: [Maintainer to add the result from the e-mail.]
* macOS builder: not run; its upload service returned HTTP 502 on
  20 November 2026.

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
