## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(mimicsurv)

## ----basic-analysis-----------------------------------------------------------
# Example data from a hypothetical study
time_points <- c(0, 6, 12, 18, 24)
n_risk <- c(200, 165, 130, 95, 65)
n_censored <- c(0, 15, 35, 53, 65)

# Extract survival analysis results
result <- extractfromKM(time_points, n_risk, n_censored)

# View hazard table
print(result$hazard_table)

# View median survival time
cat("Median survival time:", result$median_survival, "months\n")

## ----simulation---------------------------------------------------------------
# Define piecewise exponential parameters
time_points <- c(0, 6, 12, 18, 24)
hazard_rates <- c(0.05, 0.08, 0.12, 0.15)

# Simulate survival data
sim_data <- simPE(
  n = 100,  # Reduced sample size for faster execution
  time_points = time_points,
  hazard_rates = hazard_rates,
  max_time = 24,
  censoring_prob = 0.1  # Reduced censoring probability
)

# Check event rate
table(sim_data$status)

## ----summary-table------------------------------------------------------------
# Create Kaplan-Meier style summary
km_summary <- summaryKM(sim_data, time_points)

print("Number at risk:")
print(km_summary$n_risk)

print("Cumulative censored:")
print(km_summary$n_censored_cumulative)

## ----validation, eval=FALSE---------------------------------------------------
# # Run validation study (small example)
# validation_result <- validate_mimicsurv(
#   true_time_points = time_points,
#   true_hazard_rates = hazard_rates,
#   sample_sizes = c(100, 200),
#   n_simulations = 50,
#   censoring_prob = 0.15
# )
# 
# # View summary statistics
# print(validation_result$summary_statistics)

## ----true-median--------------------------------------------------------------
# Calculate true median for comparison
true_median <- getMediansurv(time_points, hazard_rates)
cat("True median survival:", round(true_median, 2), "months\n")

