# mimicsurv <img src="man/figures/mimicsurv_sticker.png" align="right" height="139" />

> Extract Survival Analysis Results from Kaplan-Meier Tables Using Person-Years Method

<!-- badges: start -->
[![R-CMD-check](https://github.com/gosukehommaEX/mimicsurv/workflows/R-CMD-check/badge.svg)](https://github.com/gosukehommaEX/mimicsurv/actions)
[![CRAN status](https://www.r-pkg.org/badges/version/mimicsurv)](https://CRAN.R-project.org/package=mimicsurv)
[![Codecov test coverage](https://codecov.io/gh/gosukehommaEX/mimicsurv/branch/main/graph/badge.svg)](https://codecov.io/gh/gosukehommaEX/mimicsurv?branch=main)
<!-- badges: end -->

## Overview

The `mimicsurv` package provides tools for survival analysis from published Kaplan-Meier tables using the person-years method. This package enables researchers to extract quantitative survival analysis results from published studies when individual patient data is not available.

It also includes simulation functions for validation of the estimation methods.

## Documentation

### Quick Start Guide

📖 **[Getting Started with mimicsurv](https://gosukehommaEX.github.io/mimicsurv/articles/getting-started.html)** - Complete tutorial with examples

### Real-World Case Studies

🏥 **[KEYNOTE-859 Clinical Trial Analysis](https://gosukehommaEX.github.io/mimicsurv/articles/keynote859-analysis.html)** - Reproducing survival results from a real pembrolizumab clinical trial

### Reference Documentation

📚 **[Package Reference](https://gosukehommaEX.github.io/mimicsurv/reference/index.html)** - Detailed function documentation

## Installation

You can install the development version of mimicsurv from GitHub:

```r
# Install devtools if you haven't already
install.packages("devtools")

# Install mimicsurv from GitHub
devtools::install_github("gosukehommaEX/mimicsurv")
```

## Dependencies

This package requires the following R packages:
- `survival` (automatically installed)
- For vignettes: `ggplot2`, `dplyr` (suggested packages)

## Features

- **Extract survival metrics** from published Kaplan-Meier curves
- **Estimate hazard rates** using the person-years method  
- **Calculate median survival times** when not explicitly reported
- **Validate methodology** with simulation functions
- **Real-world examples** including KEYNOTE-859 clinical trial analysis
- **Enable meta-analyses** and comparative effectiveness research

## Quick Example

```r
library(mimicsurv)

# Example data from a clinical trial
time_points <- c(0, 6, 12, 18, 24)
n_risk <- c(200, 150, 100, 60, 35)
n_censored <- c(0, 10, 20, 30, 40)

# Extract survival analysis results
result <- extractfromKM(time_points, n_risk, n_censored)

# View results
print(result$hazard_table)
cat("Median survival time:", result$median_survival, "months")
```

## Methodology

The package implements the person-years method for hazard estimation based on the following assumptions:

- **Piecewise exponential survival** within each time interval
- **Non-informative censoring** independent of the event process  
- **Uniform distribution** of events and censoring within intervals

For detailed mathematical background, see the [Getting Started vignette](https://gosukehommaEX.github.io/mimicsurv/articles/getting-started.html).

## Use Cases

- **Systematic reviews** and meta-analyses
- **Comparative effectiveness research**
- **Health technology assessment**
- **Clinical trial re-analysis** and validation
- **Academic research** when individual patient data (IPD) is unavailable

## Real-World Applications

The package includes a comprehensive case study analyzing the KEYNOTE-859 clinical trial, demonstrating how to:

- Extract survival data from published Kaplan-Meier curves
- Compare treatment arms (pembrolizumab vs. placebo)
- Estimate hazard rates over time
- Calculate median survival times
- Visualize results with publication-ready figures

See the [KEYNOTE-859 Analysis vignette](https://gosukehommaEX.github.io/mimicsurv/articles/keynote859-analysis.html) for the complete workflow.

## License

MIT + file LICENSE

## Citation

If you use this package in your research, please cite:

```
Homma, G. (2025). mimicsurv: Extract Survival Analysis Results from 
Kaplan-Meier Tables Using Person-Years Method. R package version 0.1.0. 
https://github.com/gosukehommaEX/mimicsurv
```
