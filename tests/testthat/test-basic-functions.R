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

# Tests for mstfromExpdists function (updated for event-based mixing)
test_that("mstfromExpdists works with basic two-group input", {
  # Basic two-group example with expected events
  expected_events <- c(75, 120)  # Expected events instead of sample sizes
  MST_subgroups <- c(6, 7.5)

  result <- mstfromExpdists(expected_events, MST_subgroups)

  # Check structure
  expect_type(result, "list")
  expect_true("subgroup_summary" %in% names(result))
  expect_true("MST_overall" %in% names(result))

  # Check subgroup_summary structure
  expect_s3_class(result$subgroup_summary, "data.frame")
  expect_equal(nrow(result$subgroup_summary), 2)
  expect_equal(ncol(result$subgroup_summary), 3)
  expect_equal(names(result$subgroup_summary), c("subgroup", "expected_events", "MST"))

  # Check values
  expect_equal(result$subgroup_summary$expected_events, expected_events)
  expect_equal(result$subgroup_summary$MST, MST_subgroups)
  expect_equal(result$subgroup_summary$subgroup, c(1, 2))

  # Check MST_overall
  expect_type(result$MST_overall, "double")
  expect_length(result$MST_overall, 1)
  expect_true(is.finite(result$MST_overall))
  expect_true(result$MST_overall > 0)

  # MST should be between the individual MSTs
  expect_true(result$MST_overall >= min(MST_subgroups))
  expect_true(result$MST_overall <= max(MST_subgroups))
})

test_that("mstfromExpdists works with three-group input", {
  # Three-group example with expected events
  expected_events <- c(150, 180, 60)  # Expected events
  MST_subgroups <- c(8, 12, 15)

  result <- mstfromExpdists(expected_events, MST_subgroups)

  # Check structure
  expect_equal(nrow(result$subgroup_summary), 3)
  expect_equal(result$subgroup_summary$expected_events, expected_events)
  expect_equal(result$subgroup_summary$MST, MST_subgroups)
  expect_equal(result$subgroup_summary$subgroup, c(1, 2, 3))

  # Check MST_overall is reasonable
  expect_true(result$MST_overall > 0)
  expect_true(result$MST_overall >= min(MST_subgroups))
  expect_true(result$MST_overall <= max(MST_subgroups))
})

test_that("mstfromExpdists handles equal expected events correctly", {
  # Equal expected events
  expected_events <- c(100, 100)
  MST_subgroups <- c(6, 12)

  result <- mstfromExpdists(expected_events, MST_subgroups)

  # With equal weights, result should be closer to arithmetic mean than extreme values
  arithmetic_mean <- mean(MST_subgroups)
  expect_true(abs(result$MST_overall - arithmetic_mean) < 2)
})

test_that("mstfromExpdists handles unequal expected events correctly", {
  # Very unequal expected events - should be closer to larger group
  expected_events <- c(10, 1000)  # 1:100 ratio
  MST_subgroups <- c(5, 20)

  result <- mstfromExpdists(expected_events, MST_subgroups)

  # Result should be much closer to the MST of the group with more events (20)
  expect_true(result$MST_overall > 15)  # Closer to 20 than to 5
})

test_that("mstfromExpdists input validation works", {
  # Test mismatched lengths
  expect_error(
    mstfromExpdists(c(100, 150), c(6)),
    "expected_events and MST_subgroups must have the same length"
  )

  # Test insufficient subgroups
  expect_error(
    mstfromExpdists(c(100), c(6)),
    "At least 2 subgroups are required"
  )

  # Test negative expected events
  expect_error(
    mstfromExpdists(c(-100, 150), c(6, 7.5)),
    "All expected event counts must be positive"
  )

  # Test zero expected events
  expect_error(
    mstfromExpdists(c(0, 150), c(6, 7.5)),
    "All expected event counts must be positive"
  )

  # Test negative MST values
  expect_error(
    mstfromExpdists(c(100, 150), c(-6, 7.5)),
    "All median survival times must be positive"
  )

  # Test zero MST values
  expect_error(
    mstfromExpdists(c(100, 150), c(0, 7.5)),
    "All median survival times must be positive"
  )

  # Test invalid search range
  expect_error(
    mstfromExpdists(c(100, 150), c(6, 7.5), search_range = c(10, 5)),
    "search_range must be a vector of length 2 with search_range\\[1\\] < search_range\\[2\\]"
  )

  # Test wrong search range length
  expect_error(
    mstfromExpdists(c(100, 150), c(6, 7.5), search_range = c(0, 50, 100)),
    "search_range must be a vector of length 2"
  )
})

# Tests for calculate_expected_events helper function
test_that("calculate_expected_events works correctly", {
  sample_sizes <- c(100, 150)
  MST_subgroups <- c(6, 7.5)
  follow_up_time <- 24

  expected_events <- calculate_expected_events(sample_sizes, MST_subgroups, follow_up_time)

  # Check output structure
  expect_type(expected_events, "double")
  expect_length(expected_events, 2)
  expect_true(all(expected_events > 0))
  expect_true(all(expected_events <= sample_sizes))  # Events can't exceed sample size

  # Check that shorter survival leads to more events
  expect_true(expected_events[1] > expected_events[2])  # MST 6 vs 7.5
})

test_that("calculate_expected_events input validation works", {
  # Test mismatched lengths
  expect_error(
    calculate_expected_events(c(100, 150), c(6)),
    "sample_sizes and MST_subgroups must have the same length"
  )

  # Test negative sample sizes
  expect_error(
    calculate_expected_events(c(-100, 150), c(6, 7.5)),
    "All sample sizes must be positive"
  )

  # Test negative MST
  expect_error(
    calculate_expected_events(c(100, 150), c(-6, 7.5)),
    "All median survival times must be positive"
  )

  # Test negative follow-up time
  expect_error(
    calculate_expected_events(c(100, 150), c(6, 7.5), -24),
    "follow_up_time must be positive"
  )
})

# Tests for mstfromExpdists_samples wrapper function
test_that("mstfromExpdists_samples wrapper works correctly", {
  sample_sizes <- c(100, 150)
  MST_subgroups <- c(6, 7.5)
  follow_up_time <- 24

  result <- mstfromExpdists_samples(sample_sizes, MST_subgroups, follow_up_time)

  # Check structure includes both sample sizes and expected events
  expect_s3_class(result$subgroup_summary, "data.frame")
  expect_true("sample_size" %in% names(result$subgroup_summary))
  expect_true("expected_events" %in% names(result$subgroup_summary))
  expect_true("follow_up_time" %in% names(result$subgroup_summary))

  # Check values
  expect_equal(result$subgroup_summary$sample_size, sample_sizes)
  expect_equal(result$subgroup_summary$MST, MST_subgroups)
  expect_true(all(result$subgroup_summary$follow_up_time == follow_up_time))

  # Check MST_overall
  expect_type(result$MST_overall, "double")
  expect_true(is.finite(result$MST_overall))
  expect_true(result$MST_overall > 0)
})

test_that("mstfromExpdists edge cases", {
  # Very similar MSTs
  expected_events <- c(100, 100)
  MST_subgroups <- c(10.0, 10.1)

  result <- mstfromExpdists(expected_events, MST_subgroups)

  # Result should be very close to both values
  expect_true(abs(result$MST_overall - 10.05) < 0.1)

  # Single dominant group by events
  expected_events_dominant <- c(1, 1000)
  MST_subgroups_dominant <- c(5, 50)

  result_dominant <- mstfromExpdists(expected_events_dominant, MST_subgroups_dominant)

  # Should be very close to the dominant group's MST
  expect_true(result_dominant$MST_overall > 45)
})

test_that("mstfromExpdists handles high MST values with appropriate search range", {
  # Very high MSTs require larger search range
  expected_events <- c(100, 100)
  MST_subgroups <- c(80, 120)  # High but reasonable MSTs

  # Should complete without error with extended search range
  result <- mstfromExpdists(expected_events, MST_subgroups, search_range = c(0, 200))
  expect_true(is.finite(result$MST_overall))
  expect_true(result$MST_overall > 0)
  expect_true(result$MST_overall >= min(MST_subgroups))
  expect_true(result$MST_overall <= max(MST_subgroups))
})

test_that("mstfromExpdists warns when median is beyond search range", {
  # Very high MSTs with default search range should warn and return NA
  expected_events <- c(100, 100)
  MST_subgroups <- c(1000, 2000)  # Very high MSTs

  # Should warn about search range and return NA
  expect_warning(
    result <- mstfromExpdists(expected_events, MST_subgroups),
    "Median survival time may be beyond the upper search range"
  )
  expect_true(is.na(result$MST_overall))
})

test_that("mstfromExpdists mathematical correctness", {
  # Test mathematical properties
  expected_events <- c(300, 200)
  MST_subgroups <- c(10, 20)

  result <- mstfromExpdists(expected_events, MST_subgroups)

  # The result should be different from simple weighted average
  simple_weighted_avg <- sum(expected_events * MST_subgroups) / sum(expected_events)
  expect_false(abs(result$MST_overall - simple_weighted_avg) < 0.001)

  # Should be between min and max
  expect_true(result$MST_overall > min(MST_subgroups))
  expect_true(result$MST_overall < max(MST_subgroups))

  # Should be closer to the group with more expected events
  # Group 1 has 300 events with MST=10, Group 2 has 200 events with MST=20
  # Result should be closer to 10 than to 20
  expect_true(abs(result$MST_overall - 10) < abs(result$MST_overall - 20))
})
