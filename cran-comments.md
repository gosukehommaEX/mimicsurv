# CRAN Submission Comments

## Test environments

* local R installation, R 4.3.2
* ubuntu 20.04 (on GitHub Actions), R 4.3.2
* win-builder (devel and release)
* macOS (on GitHub Actions), R 4.3.2

## R CMD check results

0 errors | 0 warnings | 0 notes

## Downstream dependencies

This is a new package, so there are no downstream dependencies to check.

## Additional Comments

This is the initial submission of the mimicsurv package. The package provides functions for survival analysis from Kaplan-Meier tables using the person-years method, along with simulation functions for validation.

### Package Features:
- Estimation of hazard rates from published Kaplan-Meier data
- Calculation of median survival times using piecewise exponential models
- Comprehensive simulation framework for method validation
- Integration with the PWEALL package for enhanced simulation capabilities

### Mathematical Methods:
- Person-years method for hazard rate estimation
- Piecewise exponential survival modeling
- Trapezoidal rule for person-time calculation

All functions are thoroughly documented with roxygen2, include mathematical formulations, and provide working examples.
