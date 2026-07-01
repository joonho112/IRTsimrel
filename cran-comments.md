# IRTsimrel 0.2.0 CRAN Comments

## R CMD check results

Release-candidate local check result:

- `R CMD check --as-cran --no-manual --timings`: 0 errors, 0 warnings, 1 note.
- NOTE: new submission.

Current preflight build checks:

- `R CMD build` succeeds with vignettes.
- Source tarball size: approximately 3.0 MB.
- Installed package size: approximately 5.4 MB.
- Installed `doc` directory size: approximately 4.9 MB.

## New submission

This is the first CRAN submission of `IRTsimrel`.

## R version floor

`DESCRIPTION` declares `R (>= 4.0.0)`. The GitHub Actions workflow covers R
release, devel, oldrel-1, and a separate R 4.0.5 core floor check. The floor
job disables forced suggested packages and skips vignettes because current CRAN
versions of some optional tooling packages require newer R versions. Static
compatibility checks found no native pipe (`|>`) or anonymous-function
shorthand (`\(x)`) in package R sources or tests.

## Installed size note

If CRAN reports an installed-size NOTE, the size is dominated by the package
vignettes in `inst/doc`. The package includes 12 vignettes because the main
purpose of `IRTsimrel` is methodological simulation design for
reliability-targeted item response theory studies, and the vignettes document
the two calibration algorithms, reliability conventions, applied workflows,
validation checks, and API reference material needed to use the package safely.

There are no bundled datasets or large binary assets. The largest installed
documentation files in the current preflight build are:

- `latent-distributions.html` (~1.2 MB)
- `api-reference.html` (~0.7 MB)
- `applied-guide.html` (~0.5 MB)

The documentation is retained intentionally to make the calibration assumptions,
metric choices, and validation workflow transparent to users.
