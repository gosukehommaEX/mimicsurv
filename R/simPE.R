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
#'   censoring_prob = 0.15
#' )
#'
#' # Check event rate
#' table(sim_data$status)
#'
#' @importFrom stats runif
#'
#' @export

simPE <- function(n, time_points, hazard_rates, max_time, censoring_prob = 0.1) {

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

  # Initialize output vectors
  survival_times <- numeric(n)
  event_status <- numeric(n)

  for (i in 1:n) {
    # Generate uniform random number
    u <- runif(1)

    # Find which interval the event would occur in
    survival_target <- 1 - u

    # Check if event occurs within study period
    if (survival_target <= survival_probs[length(survival_probs)]) {
      # Administrative censoring
      survival_times[i] <- max_time
      event_status[i] <- 0
    } else {
      # Find interval containing the event using manual search
      # since survival_probs is decreasing from 1 to smallest value
      interval_idx <- 1
      for (j in 1:length(survival_probs)) {
        if (survival_target > survival_probs[j]) {
          interval_idx <- j
          break
        }
      }
      interval_idx <- min(interval_idx, n_intervals)

      # Calculate exact time within interval
      if (hazard_rates[interval_idx] > 0) {
        t_start <- time_points[interval_idx]
        prev_survival <- if (interval_idx == 1) 1.0 else survival_probs[interval_idx]

        # Ensure we don't take log of negative number
        if (prev_survival > 0 && survival_target > 0 && prev_survival >= survival_target) {
          survival_time <- t_start + (1 / hazard_rates[interval_idx]) * log(prev_survival / survival_target)
        } else {
          survival_time <- max_time
        }
      } else {
        # If hazard is zero, event doesn't occur in this interval
        survival_time <- max_time
      }

      # Ensure survival time is within bounds
      survival_time <- max(0, min(survival_time, max_time))

      # Apply random censoring
      if (runif(1) < censoring_prob && survival_time < max_time) {
        # Random censoring
        censoring_time <- runif(1, 0, survival_time)
        survival_times[i] <- censoring_time
        event_status[i] <- 0
      } else {
        # Event occurs
        survival_times[i] <- survival_time
        event_status[i] <- 1
      }
    }
  }

  return(data.frame(
    time = survival_times,
    status = event_status
  ))
}
