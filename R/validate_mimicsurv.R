#' Comprehensive Simulation Study for Hazard Estimation Validation
#'
#' This function conducts a comprehensive simulation study to validate the
#' performance of hazard estimation methods using the person-years approach
#' implemented in the mimicsurv package.
#'
#' @param true_time_points Time points defining intervals for simulation
#' @param true_hazard_rates True hazard rates for each interval
#' @param sample_sizes Vector of sample sizes to test (default: c(100, 200, 500))
#' @param n_simulations Number of simulation replications (default: 100)
#' @param max_followup Maximum follow-up time (default: max of time_points)
#' @param censoring_prob Probability of random censoring (default: 0.1)
#' @param use_pweall Logical, whether to use PWEALL package if available (default: FALSE)
#'
#' @return List containing:
#'   \item{simulation_parameters}{Input parameters and true values}
#'   \item{detailed_results}{Detailed results for each sample size}
#'   \item{summary_statistics}{Summary statistics across all sample sizes}
#'
#' @details
#' The function performs the following steps for each sample size:
#' \enumerate{
#'   \item Generate simulated data from piecewise exponential distribution
#'   \item Create Kaplan-Meier style summary tables
#'   \item Apply hazard estimation method
#'   \item Calculate bias, MSE, and coverage statistics
#'   \item Summarize performance across replications
#' }
#'
#' Performance metrics include:
#' \itemize{
#'   \item Bias: Mean difference between estimated and true values
#'   \item MSE: Mean squared error
#'   \item Coverage: Proportion of estimates within acceptable range
#'   \item Convergence rate: Proportion of successful estimations
#' }
#'
#' @examples
#' # Small simulation study
#' true_times <- c(0, 6, 12, 18, 24)
#' true_hazards <- c(0.05, 0.08, 0.12, 0.15)
#'
#' validation_result <- validate_mimicsurv(
#'   true_time_points = true_times,
#'   true_hazard_rates = true_hazards,
#'   sample_sizes = c(100, 200),
#'   n_simulations = 50,
#'   censoring_prob = 0.15
#' )
#'
#' # View summary statistics
#' print(validation_result$summary_statistics)
#'
#' @importFrom stats rbinom runif
#'
#' @export

validate_mimicsurv <- function(true_time_points, true_hazard_rates,
                               sample_sizes = c(100, 200, 500),
                               n_simulations = 100,
                               max_followup = NULL,
                               censoring_prob = 0.1,
                               use_pweall = FALSE) {

  # Input validation
  if (length(true_time_points) != length(true_hazard_rates) + 1) {
    stop("Length of true_time_points should be length of true_hazard_rates + 1")
  }

  if (any(true_hazard_rates < 0)) {
    stop("All hazard rates must be non-negative")
  }

  if (is.null(max_followup)) {
    max_followup <- max(true_time_points)
  }

  # Calculate true median survival
  true_median <- getMediansurv(true_time_points, true_hazard_rates)

  # Initialize results storage
  results <- list()

  for (sample_size in sample_sizes) {
    cat(sprintf("Running simulations for sample size %d...\n", sample_size))

    # Storage for this sample size
    hazard_estimates <- array(NA, dim = c(n_simulations, length(true_hazard_rates)))
    median_estimates <- numeric(n_simulations)
    convergence_status <- logical(n_simulations)

    for (sim in 1:n_simulations) {

      tryCatch({

        # Generate data
        if (use_pweall && requireNamespace("PWEALL", quietly = TRUE)) {
          # Use PWEALL package for simulation
          sim_data_raw <- PWEALL::rpwe(
            n = sample_size,
            time = true_time_points[-1],  # Remove initial 0
            rate = true_hazard_rates
          )

          # Convert to data frame format expected by our functions
          # Apply censoring
          if (censoring_prob > 0) {
            n_to_censor <- rbinom(1, sample_size, censoring_prob)
            if (n_to_censor > 0) {
              censor_indices <- sample(sample_size, n_to_censor)
              censor_times <- runif(n_to_censor, 0, pmin(sim_data_raw, max_followup))
              sim_data_raw[censor_indices] <- censor_times
              event_status <- rep(1, sample_size)
              event_status[censor_indices] <- 0
            } else {
              event_status <- rep(1, sample_size)
            }
          } else {
            event_status <- rep(1, sample_size)
          }

          # Apply administrative censoring
          admin_censor <- sim_data_raw > max_followup
          sim_data_raw[admin_censor] <- max_followup
          event_status[admin_censor] <- 0

          sim_data <- data.frame(
            time = sim_data_raw,
            status = event_status
          )
        } else {
          # Use internal simulation function
          sim_data <- simPE(
            n = sample_size,
            time_points = true_time_points,
            hazard_rates = true_hazard_rates,
            max_time = max_followup,
            censoring_prob = censoring_prob
          )
        }

        # Create KM table
        km_table <- summaryKM(sim_data, true_time_points)

        # Apply hazard estimation method
        estimation_result <- extractfromKM(
          time_points = true_time_points,
          n_risk = km_table$n_risk,
          n_censored = km_table$n_censored_cumulative,
          warn_negative_events = FALSE
        )

        # Store results
        hazard_estimates[sim, ] <- estimation_result$hazard_table$hazard_rate
        median_estimates[sim] <- estimation_result$median_survival
        convergence_status[sim] <- TRUE

      }, error = function(e) {
        # Mark as failed convergence
        convergence_status[sim] <- FALSE
        cat(sprintf("Simulation %d failed: %s\n", sim, e$message))
      })
    }

    # Calculate performance statistics for this sample size
    valid_sims <- which(convergence_status)
    n_valid <- length(valid_sims)

    if (n_valid > 0) {

      # Hazard rate statistics
      hazard_bias <- colMeans(hazard_estimates[valid_sims, , drop = FALSE]) - true_hazard_rates
      hazard_mse <- colMeans((hazard_estimates[valid_sims, , drop = FALSE] -
                                matrix(true_hazard_rates, nrow = n_valid, ncol = length(true_hazard_rates), byrow = TRUE))^2)
      hazard_coverage <- numeric(length(true_hazard_rates))  # Placeholder for CI coverage

      # Median survival statistics
      valid_medians <- median_estimates[valid_sims]
      valid_medians <- valid_medians[!is.na(valid_medians)]

      if (length(valid_medians) > 0 && !is.na(true_median)) {
        median_bias <- mean(valid_medians) - true_median
        median_mse <- mean((valid_medians - true_median)^2)
        median_coverage_rate <- mean(abs(valid_medians - true_median) <= 0.1 * true_median)  # Within 10%
      } else {
        median_bias <- NA
        median_mse <- NA
        median_coverage_rate <- NA
      }

      # Store results for this sample size
      results[[paste0("n_", sample_size)]] <- list(
        sample_size = sample_size,
        n_valid_simulations = n_valid,
        convergence_rate = n_valid / n_simulations,

        hazard_results = data.frame(
          interval = paste0("[", true_time_points[-length(true_time_points)],
                            ",", true_time_points[-1], ")"),
          true_hazard = true_hazard_rates,
          mean_estimated_hazard = colMeans(hazard_estimates[valid_sims, , drop = FALSE]),
          bias = hazard_bias,
          mse = hazard_mse,
          stringsAsFactors = FALSE
        ),

        median_results = list(
          true_median = true_median,
          mean_estimated_median = ifelse(length(valid_medians) > 0, mean(valid_medians), NA),
          bias = median_bias,
          mse = median_mse,
          coverage_rate = median_coverage_rate,
          estimates = valid_medians
        )
      )
    }
  }

  # Overall summary
  summary_stats <- data.frame(
    sample_size = sample_sizes,
    convergence_rate = sapply(results, function(x) x$convergence_rate),
    median_bias = sapply(results, function(x) x$median_results$bias),
    median_mse = sapply(results, function(x) x$median_results$mse),
    stringsAsFactors = FALSE
  )

  return(list(
    simulation_parameters = list(
      true_time_points = true_time_points,
      true_hazard_rates = true_hazard_rates,
      true_median = true_median,
      sample_sizes = sample_sizes,
      n_simulations = n_simulations,
      max_followup = max_followup,
      censoring_prob = censoring_prob
    ),
    detailed_results = results,
    summary_statistics = summary_stats
  ))
}
