# Test basic functions in mimicsurv package

test_that("extractfromKM works with basic input", {
  # Basic test data
  time_points <- c(0, 6, 12, 18, 24)
  n_risk <- c(200, 165, 130, 95, 65)
  n_censored <- c(0, 15, 35, 53, 65)

  result <- extractfromKM(time_points, n_risk, n_censored)

  # Check structure
  expect_type(result, "list")
  expect_true("hazard_table" %in% names(result))
  expect_true("median_survival" %in% names(result))

  # Check hazard table structure
  expect_s3_class(result$hazard_table, "data.frame")
  expect_equal(nrow(result$hazard_table), 4)  # 4 intervals
  expected_cols <- c("interval", "n_at_risk_start", "n_censored_interval",
                     "n_events", "hazard_rate")
  expect_true(all(expected_cols %in% names(result$hazard_table)))

  # Check that hazard rates are non-negative
  expect_true(all(result$hazard_table$hazard_rate >= 0))

  # Check median survival
  expect_type(result$median_survival, "double")
  expect_length(result$median_survival, 1)
})

test_that("extractfromKM input validation works", {
  # Test mismatched lengths
  expect_error(
    extractfromKM(c(0, 6), c(100, 50), c(0, 10, 20)),
    "must have the same length"
  )

  # Test insufficient time points
  expect_error(
    extractfromKM(c(0), c(100), c(0)),
    "At least 2 time points are required"
  )

  # Test non-increasing time points
  expect_error(
    extractfromKM(c(0, 6, 5), c(100, 80, 60), c(0, 10, 20)),
    "must be strictly increasing"
  )

  # Test negative values
  expect_error(
    extractfromKM(c(0, 6, 12), c(-100, 80, 60), c(0, 10, 20)),
    "must be non-negative"
  )
})

test_that("getMediansurv works correctly", {
  # Test case where median is reachable
  time_points <- c(0, 6, 12, 18, 24)
  hazard_rates <- c(0.05, 0.08, 0.12, 0.15)

  median_time <- getMediansurv(time_points, hazard_rates)

  expect_type(median_time, "double")
  expect_length(median_time, 1)
  expect_true(median_time > 0)

  # Test case where median is not reachable (very low hazard)
  low_hazard_rates <- c(0.001, 0.001, 0.001, 0.001)
  median_time_na <- getMediansurv(time_points, low_hazard_rates)
  expect_true(is.na(median_time_na))
})

test_that("getMediansurv input validation works", {
  # Test mismatched lengths
  expect_error(
    getMediansurv(c(0, 6, 12), c(0.05, 0.08, 0.12, 0.15)),
    "Length of time_points should be length of hazard_rates \\+ 1"
  )

  # Test negative hazard rates
  expect_error(
    getMediansurv(c(0, 6, 12), c(-0.05, 0.08)),
    "All hazard rates must be non-negative"
  )
})

test_that("simPE generates appropriate data", {
  time_points <- c(0, 6, 12, 18, 24)
  hazard_rates <- c(0.05, 0.08, 0.12, 0.15)

  sim_data <- simPE(
    n = 100,
    time_points = time_points,
    hazard_rates = hazard_rates,
    max_time = 24,
    censoring_prob = 0.1
  )

  # Check structure
  expect_s3_class(sim_data, "data.frame")
  expect_equal(nrow(sim_data), 100)
  expect_true(all(c("time", "status") %in% names(sim_data)))

  # Check data ranges
  expect_true(all(sim_data$time >= 0))
  expect_true(all(sim_data$time <= 24))
  expect_true(all(sim_data$status %in% c(0, 1)))

  # Should have some censored observations with censoring_prob > 0
  expect_true(any(sim_data$status == 0))
})

test_that("summaryKM creates proper summary", {
  # Create simple test data
  sim_data <- data.frame(
    time = c(5, 10, 15, 20, 25),
    status = c(1, 0, 1, 1, 0)
  )
  time_points <- c(0, 6, 12, 18, 24)

  km_summary <- summaryKM(sim_data, time_points)

  # Check structure
  expect_type(km_summary, "list")
  expect_true("n_risk" %in% names(km_summary))
  expect_true("n_censored_cumulative" %in% names(km_summary))

  # Check lengths
  expect_equal(length(km_summary$n_risk), length(time_points))
  expect_equal(length(km_summary$n_censored_cumulative), length(time_points))

  # Check that n_risk is decreasing
  expect_true(all(diff(km_summary$n_risk) <= 0))

  # Check that cumulative censoring is increasing
  expect_true(all(diff(km_summary$n_censored_cumulative) >= 0))
})

test_that("summaryKM input validation works", {
  # Test non-data.frame input
  expect_error(
    summaryKM(list(time = c(1, 2), status = c(1, 0)), c(0, 5)),
    "sim_data must be a data frame"
  )

  # Test missing columns
  expect_error(
    summaryKM(data.frame(x = c(1, 2), y = c(1, 0)), c(0, 5)),
    "must contain 'time' and 'status' columns"
  )
})
