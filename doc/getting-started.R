## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = TRUE,
  echo = TRUE
)

## ----setup--------------------------------------------------------------------
library(mimicsurv)

## ----basic-analysis-----------------------------------------------------------
# Example data from a hypothetical study
time_points <- c(0, 6, 12, 18, 24)
n_risk <- c(200, 150, 100, 60, 35)
n_censored <- c(0, 10, 20, 30, 40)

# Extract survival analysis results
result <- extractfromKM(time_points, n_risk, n_censored)

# View hazard table
print(result$hazard_table)

# View median survival time
cat("Median survival time:", result$median_survival, "months\n")

## ----mixed-exponential-basic--------------------------------------------------
# Example: Two biomarker subgroups with known expected events
# In practice, these would come from published studies or trial design
expected_events <- c(75, 120)  # Expected events from each subgroup
MST_subgroups <- c(6, 7.5)     # Median survival times from each subgroup

print(paste("Expected events:", paste(expected_events, collapse = ", ")))

# Use expected events for mixing probabilities
result_mixed <- mstfromExpdists(expected_events, MST_subgroups)

# View subgroup summary
print(result_mixed$subgroup_summary)

# View overall median survival time
cat("Overall MST:", round(result_mixed$MST_overall, 2), "months\n")

## ----mixed-exponential-three--------------------------------------------------
# Example: Three treatment subgroups with known expected events
expected_events_3 <- c(150, 180, 60)  # Expected events from published data
MST_subgroups_3 <- c(8, 12, 15)       # Median survival times from each subgroup

result_3groups <- mstfromExpdists(expected_events_3, MST_subgroups_3)

# View results
print(result_3groups$subgroup_summary)
cat("Overall MST:", round(result_3groups$MST_overall, 2), "months\n")

# Compare with simple weighted average (incorrect method)
# For comparison: if we incorrectly used sample sizes
sample_sizes_3 <- c(200, 300, 100)  # Hypothetical sample sizes
weighted_avg <- sum(sample_sizes_3 * MST_subgroups_3) / sum(sample_sizes_3)
cat("Simple weighted average:", round(weighted_avg, 2), "months\n")
cat("Difference:", round(abs(result_3groups$MST_overall - weighted_avg), 2), "months\n")

## ----mixed-exponential-wrapper------------------------------------------------
# When only sample sizes are available from published studies
sample_sizes <- c(100, 150)
MST_subgroups <- c(6, 7.5)
follow_up_time <- 24  # months

# Calculate expected events using exponential distribution formula
lambda <- log(2) / MST_subgroups
expected_events <- sample_sizes * (1 - exp(-lambda * follow_up_time))

result_from_samples <- mstfromExpdists(expected_events, MST_subgroups)

print(data.frame(
  subgroup = result_from_samples$subgroup_summary$subgroup,
  sample_size = sample_sizes,
  expected_events = round(expected_events, 1),
  MST = MST_subgroups
))
cat("Overall MST (from sample sizes):", round(result_from_samples$MST_overall, 2), "months\n")

## ----simulation---------------------------------------------------------------
# Define true parameters
true_time_points <- c(0, 6, 12, 18, 24)
true_hazard_rates <- c(0.05, 0.08, 0.12, 0.15)

# Generate simulated data
set.seed(123)
sim_data <- simPE(
  n = 100,  # Reduced sample size for faster execution
  time_points = true_time_points,
  hazard_rates = true_hazard_rates,
  max_time = 24,
  censoring_prob = 0.1,  # Reduced censoring probability
  seed = 123  # For reproducible results
)

# Check event rate
print("Event status distribution:")
table(sim_data$status)

# View first few rows
print("First 10 observations:")
head(sim_data, 10)

## ----summary-tables-----------------------------------------------------------
# Create Kaplan-Meier style summary
km_summary <- summaryKM(sim_data, true_time_points)

print("Number at risk at each time point:")
print(km_summary$n_risk)

print("Cumulative censoring at each time point:")
print(km_summary$n_censored_cumulative)

# Create comparison table
comparison_table <- data.frame(
  Time = true_time_points,
  True_n_risk = c(100, 75, 50, 25, 10),  # Hypothetical true values
  Simulated_n_risk = km_summary$n_risk,
  Difference = km_summary$n_risk - c(100, 75, 50, 25, 10)
)

print("Comparison of simulated vs. expected values:")
print(comparison_table)

## ----reproducibility----------------------------------------------------------
# Demonstrate reproducibility with seeds
sim1 <- simPE(n = 10, time_points = c(0, 12, 24), 
              hazard_rates = c(0.1, 0.2), max_time = 24, seed = 456)
sim2 <- simPE(n = 10, time_points = c(0, 12, 24), 
              hazard_rates = c(0.1, 0.2), max_time = 24, seed = 456)
sim3 <- simPE(n = 10, time_points = c(0, 12, 24), 
              hazard_rates = c(0.1, 0.2), max_time = 24, seed = 789)

print("First simulation (first 5 rows):")
print(head(sim1, 5))

print("Second simulation (first 5 rows):")
print(head(sim2, 5))

print(paste("Are the simulations identical?", identical(sim1, sim2)))
print(paste("Are sim1 and sim3 different?", !identical(sim1, sim3)))

## ----method-comparison--------------------------------------------------------
# Apply our method to simulated data
sim_result <- extractfromKM(
  time_points = true_time_points,
  n_risk = km_summary$n_risk,
  n_censored = km_summary$n_censored_cumulative
)

print("Analysis of simulated data:")
print(sim_result$hazard_table)

# Compare estimated vs. true hazard rates
comparison_hazards <- data.frame(
  Interval = sim_result$hazard_table$interval,
  True_Hazard = true_hazard_rates,
  Estimated_Hazard = sim_result$hazard_table$hazard_rate,
  Relative_Error = abs(sim_result$hazard_table$hazard_rate - true_hazard_rates) / true_hazard_rates * 100
)

print("Hazard rate comparison:")
print(comparison_hazards)

# Calculate true median survival time
true_median <- getMediansurv(true_time_points, true_hazard_rates)
print(paste("True median survival:", round(true_median, 2), "months"))
print(paste("Estimated median survival:", round(sim_result$median_survival, 2), "months"))

## ----published-data-----------------------------------------------------------
# Compare our method with published Kaplan-Meier data
published_times <- c(0, 3, 6, 9, 12, 15, 18, 21, 24)
published_n_risk <- c(180, 165, 145, 128, 110, 95, 78, 65, 50)
published_censored <- c(0, 5, 15, 22, 35, 48, 62, 75, 90)

# Apply our method
published_result <- extractfromKM(published_times, published_n_risk, published_censored)

print("Analysis of published data:")
print(published_result$hazard_table)
print(paste("Median survival:", round(published_result$median_survival, 2), "months"))

## ----efficiency---------------------------------------------------------------
# Benchmark computational performance
system.time({
  result_timing <- extractfromKM(time_points, n_risk, n_censored)
})

# Test with larger datasets
large_times <- seq(0, 60, by = 3)
large_n_risk <- seq(1000, 100, length.out = length(large_times))
large_censored <- cumsum(c(0, rep(20, length(large_times) - 1)))

print("Performance with larger dataset:")
system.time({
  large_result <- extractfromKM(large_times, large_n_risk, large_censored)
})

print("Large dataset median survival:")
print(paste(round(large_result$median_survival, 2), "months"))

## ----advanced-median----------------------------------------------------------
# Example with complex piecewise exponential parameters
complex_times <- c(0, 3, 6, 9, 12, 15, 18, 21, 24)
complex_hazards <- c(0.02, 0.05, 0.08, 0.12, 0.15, 0.18, 0.22, 0.25)

complex_median <- getMediansurv(complex_times, complex_hazards)
print(paste("Complex scenario median:", round(complex_median, 2), "months"))

## ----sensitivity--------------------------------------------------------------
# Sensitivity to censoring assumptions
censoring_levels <- c(0.05, 0.10, 0.15, 0.20)
median_estimates <- numeric(length(censoring_levels))

for (i in seq_along(censoring_levels)) {
  sim_sens <- simPE(n = 200, time_points = true_time_points, 
                    hazard_rates = true_hazard_rates, max_time = 24,
                    censoring_prob = censoring_levels[i], seed = 100 + i)
  km_sens <- summaryKM(sim_sens, true_time_points)
  result_sens <- extractfromKM(true_time_points, km_sens$n_risk, 
                               km_sens$n_censored_cumulative)
  median_estimates[i] <- result_sens$median_survival
}

sensitivity_table <- data.frame(
  Censoring = paste0(censoring_levels * 100, "%"),
  Median_Estimate = round(median_estimates, 2),
  Bias = round(median_estimates - true_median, 2)
)

print("Sensitivity to censoring rate:")
print(sensitivity_table)

