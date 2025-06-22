# mimicsurv 0.1.0

## Initial release

### New Features

* Added `extractfromKM()` function for extracting survival analysis results from Kaplan-Meier tables using person-years method
* Added `getMediansurv()` function for median survival time estimation
* Added simulation functions (`simPE()`, `summaryKM()`, `validate_mimicsurv()`) for methodology validation
* Comprehensive vignettes including:
  - Getting Started guide with theoretical foundations
  - KEYNOTE-859 clinical trial case study
* Full pkgdown website with documentation and examples

### Documentation

* Complete function reference documentation
* Mathematical background and methodology explanation
* Real-world clinical trial analysis example
* GitHub Pages website: https://gosukehommaEX.github.io/mimicsurv/

### Dependencies

* Suggests: ggplot2, dplyr for visualization in vignettes
* Suggests: PWEALL for simulation functions
* Base R implementation for core functionality
