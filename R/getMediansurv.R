#' Calculate Median Survival Time from Piecewise Exponential Parameters
#'
#' This function calculates median survival time from time points and
#' hazard rates using piecewise exponential survival function.
#'
#' @param time_points Numeric vector of time points defining intervals
#' @param hazard_rates Numeric vector of hazard rates for each interval
#'
#' @return Numeric scalar of median survival time (NA if not reached)
#'
#' @details
#' The function assumes piecewise exponential survival with constant hazard
#' within each interval. It calculates survival probabilities at each time
#' point and finds the time where survival probability equals 0.5.
#'
#' Mathematical formulation:
#' \deqn{S(t) = \exp\left(-\sum_{j=1}^{k} \lambda_j \cdot \Delta t_j\right)}
#' where \eqn{\lambda_j} is hazard rate in interval j and \eqn{\Delta t_j}
#' is the length of interval j up to time t.
#'
#' @examples
#' # Calculate median for 4-interval model
#' time_points <- c(0, 6, 12, 18, 24)
#' hazard_rates <- c(0.05, 0.08, 0.12, 0.15)
#'
#' median_time <- getMediansurv(time_points, hazard_rates)
#' print(paste("Median survival:", round(median_time, 2), "months"))
#'
#' @export

getMediansurv <- function(time_points, hazard_rates) {

  # Input validation
  if (length(time_points) != length(hazard_rates) + 1) {
    stop("Length of time_points should be length of hazard_rates + 1")
  }

  if (any(hazard_rates < 0)) {
    stop("All hazard rates must be non-negative")
  }

  n_intervals <- length(time_points) - 1
  interval_lengths <- diff(time_points)

  # Calculate survival probabilities at each time point
  survival_probs <- numeric(length(time_points))
  survival_probs[1] <- 1.0

  cumulative_hazard <- 0
  for (i in 1:n_intervals) {
    cumulative_hazard <- cumulative_hazard + hazard_rates[i] * interval_lengths[i]
    survival_probs[i + 1] <- exp(-cumulative_hazard)
  }

  # Check if median is reachable
  final_survival <- survival_probs[length(survival_probs)]
  if (final_survival >= 0.5) {
    return(NA)  # Median not reached
  }

  # Find interval containing median
  median_interval <- which(survival_probs < 0.5)[1] - 1
  if (median_interval == 0) median_interval <- 1

  # Calculate median within the interval
  t_start <- time_points[median_interval]
  lambda <- hazard_rates[median_interval]
  prev_survival <- if (median_interval == 1) 1.0 else survival_probs[median_interval]

  if (lambda == 0) {
    return(NA)  # Cannot calculate if hazard is zero
  }

  median_time <- t_start + (1 / lambda) * log(prev_survival / 0.5)

  return(median_time)
}
