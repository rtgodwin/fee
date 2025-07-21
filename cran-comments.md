This is a new package submission to CRAN. The package `fee` provides methods to estimate the first-exposure effect (FEE) using a one-inflated positive Poisson model or a one-inflated zero-truncated negative binomial model. It also estimates marginal FEEs and standard errors for both the FEE and marginal FEEs.

## CRAN checks

All checks were performed using the following:

- `devtools::check()` on local machine — **no errors, warnings, or notes**
- `devtools::check_win_devel()` — **no errors or warnings; 1 NOTE**:

Maintainer: 'Ryan T. Godwin ryan.godwin@umanitoba.ca'
New submission

This NOTE is expected for new submissions.

- `rhub::check()` on multiple platforms via GitHub Actions:

✔ Linux (ubuntu-latest)  
✔ Windows (windows-latest)  
✔ macOS (macos-13, macos-15, and macos-arm64)

All RHUB checks completed successfully with **no errors, warnings, or notes**.

## Additional notes

- I have verified that the package passes all checks with R-devel.
- The package does not require compilation.
- All dependencies are available on CRAN.