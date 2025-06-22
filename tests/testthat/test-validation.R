# Test validation functions in mimicsurv package

test_that("validate_mimicsurv works with basic input", {
  # Skip if PWEALL is not available
  skip_if_not_installed("PWEALL")

  # Small test case
  true_times <- c(0, 6, 12)
  true_hazards <- c(0.1, 0.15)

  # Run very small validation
  result <- validate_mimicsurv(
    true_time_points = true_times,
    true_hazard_rates = true_hazards,
    sample_sizes = c(50),
    n_simulations = 5,
    max_followup = 12,
    censoring_prob = 0.1,
    use_pweall = FALSE  # Use internal simulation
  )

  # Check structure
  expect_type(result, "list")
  expect_true("simulation_parameters" %in% names(result))
  expect_true("detailed_results" %in% names(result))
  expect_true("summary_statistics" %in% names(result))

  # Check simulation parameters
  params <- result$simulation_parameters
  expect_equal(params$true_time_points, true_times)
  expect_equal(params$true_hazard_rates, true_hazards)
  expect_equal(params$sample_sizes, c(50))
  expect_equal(params$n_simulations, 5)

  # Check summary statistics structure
  summary_stats <- result$summary_statistics
  expect_s3_class(summary_stats, "data.frame")
  expect_true("sample_size" %in% names(summary_stats))
  expect_true("convergence_rate" %in% names(summary_stats))
})

test_that("validate_mimicsurv input validation works", {
  # Test mismatched lengths
  expect_error(
    validate_mimicsurv(
      true_time_points = c(0, 6, 12),
      true_hazard_rates = c(0.1, 0.15, 0.2),  # Too many hazard rates
      sample_sizes = c(50),
      n_simulations = 5
    ),
    "Length of true_time_points should be length of true_hazard_rates \\+ 1"
  )

  # Test negative hazard rates
  expect_error(
    validate_mimicsurv(
      true_time_points = c(0, 6, 12),
      true_hazard_rates = c(-0.1, 0.15),
      sample_sizes = c(50),
      n_simulations = 5
    ),
    "All hazard rates must be non-negative"
  )
})

test_that("validate_mimicsurv handles PWEALL integration", {
  # Skip if PWEALL is not available
  skip_if_not_installed("PWEALL")

  true_times <- c(0, 6, 12)
  true_hazards <- c(0.1, 0.15)

  # Test with PWEALL enabled - use very small parameters to avoid issues
  result_pweall <- validate_mimicsurv(
    true_time_points = true_times,
    true_hazard_rates = true_hazards,
    sample_sizes = c(20),  # Very small sample size
    n_simulations = 2,     # Very few simulations
    max_followup = 12,
    censoring_prob = 0.05, # Low censoring
    use_pweall = TRUE
  )

  # Should complete without error
  expect_type(result_pweall, "list")
  expect_true("simulation_parameters" %in% names(result_pweall))
})

test_that("validate_mimicsurv handles missing PWEALL gracefully", {
  true_times <- c(0, 6, 12)
  true_hazards <- c(0.1, 0.15)

  # Test with PWEALL disabled (should always work)
  result_no_pweall <- validate_mimicsurv(
    true_time_points = true_times,
    true_hazard_rates = true_hazards,
    sample_sizes = c(30),
    n_simulations = 3,
    max_followup = 12,
    censoring_prob = 0.1,
    use_pweall = FALSE
  )

  # Should complete without error
  expect_type(result_no_pweall, "list")
  expect_true("simulation_parameters" %in% names(result_no_pweall))
})
