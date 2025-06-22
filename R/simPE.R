#' Simulate Survival Data from Piecewise Exponential Distribution
#'
#' This function generates simulated survival data following a piecewise
#' exponential distribution with specified time intervals and hazard rates.
#'
#' @param n Sample size (number of subjects to simulate)
#' @param time_points Numeric vector of time points defining intervals
#' @param hazard_rates Numeric vector of hazard rates for each interval
#' @param max_time Maximum follow-up time (administrative censoring)
#' @param censoring_prob Probability of random censoring (default 0.1)
#' @param seed Random seed for reproducibility (default NULL)
#'
#' @return Data frame with columns:
#'   \item{time}{Observed time (either event time or censoring time)}
#'   \item{status}{Event indicator (1 = event, 0 = censored)}
#'
#' @details
#' The function simulates survival times using inverse transform sampling
#' for piecewise exponential distribution. Random censoring is applied
#' uniformly between 0 and the simulated survival time.
#'
#' @examples
#' # Simulate data with 4 intervals
#' time_points <- c(0, 6, 12, 18, 24)
#' hazard_rates <- c(0.05, 0.08, 0.12, 0.15)
#'
#' sim_data <- simPE(
#'   n = 200,
#'   time_points = time_points,
#'   hazard_rates = hazard_rates,
#'   max_time = 24,
#'   censoring_prob = 0.15,
#'   seed = 123
#' )
#'
#' # Check event rate
#' table(sim_data$status)
#'
#' # Reproducible results
#' sim_data1 <- simPE(n = 100, time_points, hazard_rates, 24, seed = 456)
#' sim_data2 <- simPE(n = 100, time_points, hazard_rates, 24, seed = 456)
#' identical(sim_data1, sim_data2)  # Should be TRUE
#'
#' @importFrom stats runif
#'
#' @export

simPE <- function(n, time_points, hazard_rates, max_time, censoring_prob = 0.1, seed = NULL) {

  # Set seed for reproducibility if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Input validation
  if (length(time_points) != length(hazard_rates) + 1) {
    stop("Length of time_points should be length of hazard_rates + 1")
  }

  if (any(hazard_rates < 0)) {
    stop("All hazard rates must be non-negative")
  }

  n_intervals <- length(hazard_rates)
  interval_lengths <- diff(time_points)

  # Pre-calculate cumulative hazards and survival probabilities at interval boundaries
  cumulative_hazards <- numeric(length(time_points))
  survival_probs <- numeric(length(time_points))
  survival_probs[1] <- 1.0

  for (i in 1:n_intervals) {
    cumulative_hazards[i + 1] <- cumulative_hazards[i] + hazard_rates[i] * interval_lengths[i]
    survival_probs[i + 1] <- exp(-cumulative_hazards[i + 1])
  }

  # Generate survival times using inverse transform sampling
  survival_times <- numeric(n)

  for (j in 1:n) {
    u <- runif(1)
    target_cum_hazard <- -log(u)

    # Find which interval the target cumulative hazard falls into
    survival_time <- max_time  # Default to max_time if beyond all intervals

    for (i in 1:n_intervals) {
      if (target_cum_hazard <= cumulative_hazards[i + 1]) {
        # Target falls in interval i
        remaining_hazard <- target_cum_hazard - cumulative_hazards[i]
        time_in_interval <- remaining_hazard / hazard_rates[i]
        survival_time <- time_points[i] + time_in_interval
        break
      }
    }

    # If target_cum_hazard is beyond the last interval, extrapolate
    if (target_cum_hazard > cumulative_hazards[n_intervals + 1]) {
      remaining_hazard <- target_cum_hazard - cumulative_hazards[n_intervals + 1]
      time_beyond_last <- remaining_hazard / hazard_rates[n_intervals]
      survival_time <- time_points[n_intervals + 1] + time_beyond_last
    }

    survival_times[j] <- survival_time
  }

  # Apply administrative censoring
  survival_times <- pmin(survival_times, max_time)

  # Apply random censoring
  censoring_times <- runif(n, 0, max_time)
  should_censor <- runif(n) < censoring_prob

  observed_times <- ifelse(should_censor & censoring_times < survival_times,
                           censoring_times, survival_times)
  status <- ifelse(should_censor & censoring_times < survival_times, 0, 1)

  return(data.frame(time = observed_times, status = status))
}
