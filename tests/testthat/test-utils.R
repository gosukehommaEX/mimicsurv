test_that("getMedianSurv works correctly", {
  time_points <- c(0, 6, 12, 18, 24)
  hazard_rates <- c(0.05, 0.08, 0.12, 0.15)

  median_surv <- getMedianSurv(time_points, hazard_rates)

  # Check result is numeric
  expect_type(median_surv, "double")

  # Check that survival at median time is 0.5 (approximately)
  survFun <- function(t) {
    total_hazard <- 0
    for (i in seq_along(hazard_rates)) {
      if (t <= time_points[i]) {
        break
      }
      interval_length <- min(t, time_points[i+1]) - time_points[i]
      total_hazard <- total_hazard + hazard_rates[i] * interval_length
    }
    return(exp(-total_hazard))
  }

  expect_equal(survFun(median_surv), 0.5, tolerance = 1e-6)
})

test_that("simPE generates data correctly", {
  time_points <- c(0, 6, 12, 18, 24)
  hazard_rates <- c(0.05, 0.08, 0.12, 0.15)
  n <- 100
  max_time <- 36

  set.seed(123)  # For reproducibility
  data <- simPE(n, time_points, hazard_rates, max_time, censoring_prob = 0)

  # Check data structure
  expect_equal(nrow(data), n)
  expect_equal(ncol(data), 3)
  expect_true(all(c("id", "time", "status") %in% names(data)))

  # Check data types
  expect_type(data$id, "integer")
  expect_type(data$time, "double")
  expect_type(data$status, "double")

  # Check ranges
  expect_true(all(data$time >= 0))
  expect_true(all(data$status %in% c(0, 1)))

  # Check censoring
  expect_true(all(data$time[data$status == 0] <= max_time))
})

test_that("summaryKM creates correct summary table", {
  time_points <- c(0, 6, 12, 18, 24)
  hazard_rates <- c(0.05, 0.08, 0.12, 0.15)
  n <- 100

  set.seed(123)  # For reproducibility
  data <- simPE(n, time_points, hazard_rates, max_time = 36, censoring_prob = 0)

  km_table <- summaryKM(data, time_points)

  # Check structure
  expect_type(km_table, "list")
  expect_true(all(c("n_risk", "n_event", "n_censored", "n_censored_cumulative") %in% names(km_table)))

  # Check dimensions
  expect_length(km_table$n_risk, length(time_points))
  expect_length(km_table$n_event, length(time_points) - 1)
  expect_length(km_table$n_censored, length(time_points) - 1)
  expect_length(km_table$n_censored_cumulative, length(time_points))

  # Check n_risk is decreasing
  expect_true(all(diff(km_table$n_risk) <= 0))

  # Check first n_risk equals total sample
  expect_equal(km_table$n_risk[1], n)
})

test_that("extractFromKM estimates hazard rates correctly", {
  time_points <- c(0, 6, 12, 18, 24)
  hazard_rates <- c(0.05, 0.08, 0.12, 0.15)
  n <- 1000  # Large sample to reduce variation

  set.seed(123)  # For reproducibility
  data <- simPE(n, time_points, hazard_rates, max_time = 36, censoring_prob = 0)

  km_table <- summaryKM(data, time_points)

  result <- extractFromKM(
    time_points = time_points,
    n_risk = km_table$n_risk,
    n_censored = km_table$n_censored_cumulative
  )

  # Check structure
  expect_type(result, "list")
  expect_true(all(c("hazard_table", "median_survival", "survival_function") %in% names(result)))

  # Check hazard_table
  expect_s3_class(result$hazard_table, "data.frame")
  expect_equal(nrow(result$hazard_table), length(time_points) - 1)

  # Check hazard rates (should be close to true values with large sample)
  estimated_hazards <- result$hazard_table$hazard_rate

  # Allow some variation due to randomness
  for (i in seq_along(hazard_rates)) {
    expect_equal(estimated_hazards[i], hazard_rates[i], tolerance = 0.2)
  }
})
