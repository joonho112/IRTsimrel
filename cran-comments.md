# IRTsimrel 0.3.0 CRAN Comments

## Test environments

- Local macOS arm64, R 4.6.0: full `R CMD check --as-cran --timings`,
  including manual, tests, examples, and rebuilding all 12 vignettes.
- Local macOS arm64, R 4.5.1: the exact source tarball checked with
  `_R_CHECK_FORCE_SUGGESTS_=false`, `--no-manual`, and
  `--ignore-vignettes` after the full R 4.6.0 build/check.
- GitHub Actions is configured for macOS and Windows R-release, Ubuntu
  R-devel/release/oldrel-1, and a separate Ubuntu R 4.0.5 core-floor job.
  Those remote jobs must be green on the final release commit before
  submission; their results are not represented as completed local checks.

## R CMD check results

Local R 4.6.0 full check:

- 0 errors, 0 warnings, 2 notes.
- NOTE 1: new submission.
- NOTE 2: formal HTML-manual validation was skipped because the system
  `/usr/bin/tidy` is the Apple build from 2006. The HTML and PDF manuals were
  built, and package installation, examples (including `--run-donttest`),
  tests, URLs, and all vignette builds completed successfully.

Local R 4.5.1 check of the same tarball:

- 0 errors, 0 warnings, 0 notes (`Status: OK`).
- 1,746 test expectations passed; one environment-conditional IRW-absence
  test was skipped because `irw` was installed. Its absence/error path was
  also exercised separately with an isolated namespace mock.

The checked source tarball has SHA-256
`acf92436c88338d05f34361179211b0f751056ee0dca2fbc77f61bbff3e5a48e`,
size 3,279,580 bytes, and 135 archive entries.

## New submission

This is the first CRAN submission of `IRTsimrel`.

## Changes in 0.3.0

- Added end-to-end unidimensional 3PL support with logistic scaling `D = 1`.
  Item guessing is named `guessing`; `c` remains the global discrimination
  multiplier and never rescales guessing.
- Added stable extreme-predictor probability/information kernels and
  guessing-aware item generation, EQC, SAC, response simulation, and S3
  methods.
- Replaced global-monotonicity assumptions with an adaptive log-scale topology
  scan and explicit root-selection policies.
- Added SAC topology preflight, branch guards, independent evaluation, and
  explicit fixed-form versus item-superpopulation estimands.
- Preserved the historical SAC argument prefix and deprecated SPC aliases;
  legacy Rasch/2PL results without guessing are interpreted lazily as
  zero-guessing designs.

The release remains unidimensional. The public TAM helper continues to support
its documented Rasch/2PL WLE and EAP paths; external 3PL validation was EAP
only, and no TAM 3PL WLE support is claimed.

## R version floor

`DESCRIPTION` declares `R (>= 4.0.0)`. The workflow includes an R 4.0.5 core
job with forced suggested packages disabled and vignette building omitted,
because current optional documentation packages may require newer R versions.
Static compatibility checks found no native pipe (`|>`) or anonymous-function
shorthand (`\(x)`) in package R sources, tests, or vignette code. The remote
R 4.0.5 job remains a release-time gate.

## Installed size

The local full check reported an installed size of 5.8 MB, including 5.0 MB in
`doc`. The package intentionally includes 12 methodological vignettes covering
the model/estimand contracts, two calibration algorithms, applied workflows,
validation evidence, and API reference. No bundled datasets or large binary
research artifacts are shipped; development evidence, logs, and the pkgdown
site are excluded from the source tarball.

## Downstream dependencies

There are currently no known reverse dependencies.
