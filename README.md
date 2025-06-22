# mimicsurv

An R package for survival analysis from Kaplan-Meier tables using the person-years method.

## Overview

The `mimicsurv` package provides functions to estimate hazard rates and median survival time from published Kaplan-Meier table data using the person-years method, assuming exponential distribution within each interval. It also includes simulation functions for validation of the estimation methods.

## Documentation

### Quick Start Guide
📖 **[Getting Started with mimicsurv](https://gosukehommaEX.github.io/mimicsurv/articles/getting-started.html)** - Complete tutorial with examples

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
- `stats` (base R)
- `utils` (base R)
- `PWEALL` (for simulation validation)

## Main Functions

### Core Analysis Functions

- `extractfromKM()`: Extract survival analysis results from Kaplan-Meier table using person-years method
- `getMediansurv()`: Calculate median survival time from piecewise exponential parameters

### Simulation Functions

- `simPE()`: Simulate survival data from piecewise exponential distribution
- `summaryKM()`: Create Kaplan-Meier style summary table from simulated data
- `validate_mimicsurv()`: Comprehensive simulation study for hazard estimation validation

## Example Usage

```r
library(mimicsurv)

# Example: Basic survival analysis from KM table
time_points <- c(0, 6, 12, 18, 24)
n_risk <- c(200, 150, 100, 60, 35)
n_censored <- c(0, 10, 20, 30, 40)

result <- extractfromKM(time_points, n_risk, n_censored)
print(result$hazard_table)
print(paste("Median survival:", result$median_survival, "months"))

# Example: Simulation validation
true_times <- c(0, 6, 12, 18, 24)
true_hazards <- c(0.05, 0.08, 0.12, 0.15)

validation_result <- validate_mimicsurv(
  true_time_points = true_times,
  true_hazard_rates = true_hazards,
  sample_sizes = c(100, 200),
  n_simulations = 50,
  censoring_prob = 0.15
)

print(validation_result$summary_statistics)
```

## Mathematical Background

The package implements survival analysis methods based on:

- **Person-years method**: Calculates hazard rates as events per person-time using trapezoidal rule
- **Piecewise exponential model**: Assumes constant hazard within intervals
- **Median survival estimation**: Uses exponential survival function: S(t) = exp(-∑λⱼΔtⱼ)

## Author

**Gosuke Homma**  
GitHub: [gosukehommaEX](https://github.com/gosukehommaEX)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
