#' Extract Survival Analysis Results from Kaplan-Meier Table Using Person-Years Method
#'
#' This function estimates hazard rates and median survival time from published
#' Kaplan-Meier table data using the person-years method, assuming exponential
#' distribution within each interval.
#'
#' @param time_points Numeric vector of time points (length n+1 for n intervals)
#' @param n_risk Numeric vector of number at risk at each time point (length n+1)
#' @param n_censored Numeric vector of cumulative number censored at each time point (length n+1)
#' @param warn_negative_events Logical, whether to warn about negative events (default TRUE)
#'
#' @return List containing:
#'   \item{hazard_table}{Data frame with interval summaries and hazard rates}
#'   \item{median_survival}{Estimated median survival time (NA if not reached)}
#'
#' @details
#' The function calculates:
#' \itemize{
#'   \item Number of events in each interval (accounting for risk set changes)
#'   \item Person-time using trapezoidal rule
#'   \item Hazard rates as events/person-time
#'   \item Median survival time using exponential survival function
#' }
#'
#' @examples
#' # Example 1: Basic usage
#' time_points <- c(0, 6, 12, 18, 24)
#' n_risk <- c(200, 165, 130, 95, 65)
#' n_censored <- c(0, 15, 35, 53, 65)
#'
#' result <- extractfromKM(time_points, n_risk, n_censored)
#' print(result$hazard_table)
#' print(paste("Median survival:", result$median_survival, "months"))
#'
#' @importFrom stats diff
#'
#' @export

extractfromKM <- function(time_points, n_risk, n_censored,
                          warn_negative_events = TRUE) {

  # Input validation
  if (length(time_points) != length(n_risk) || length(time_points) != length(n_censored)) {
    stop("time_points, n_risk, and n_censored must have the same length")
  }

  if (length(time_points) < 2) {
    stop("At least 2 time points are required")
  }

  if (any(diff(time_points) <= 0)) {
    stop("time_points must be strictly increasing")
  }

  if (any(n_risk < 0) || any(n_censored < 0)) {
    stop("n_risk and n_censored must be non-negative")
  }

  n_intervals <- length(time_points) - 1

  # Create interval labels
  interval_labels <- paste0("[", time_points[1:n_intervals], ",",
                            time_points[2:(n_intervals + 1)], ")")

  # Calculate number censored in each interval
  n_censored_interval <- diff(n_censored)

  # Calculate number of events in each interval
  n_events <- numeric(n_intervals)
  for (i in 1:n_intervals) {
    n_events[i] <- n_risk[i] - n_risk[i + 1] - n_censored_interval[i]
  }

  # Handle negative events
  if (any(n_events < 0)) {
    if (warn_negative_events) {
      warning("Negative events detected in some intervals. Setting to 0.")
    }
    n_events[n_events < 0] <- 0
  }

  # Calculate person-time using trapezoidal rule
  person_time <- (n_risk[1:n_intervals] + n_risk[2:(n_intervals + 1)]) / 2 * diff(time_points)

  # Calculate hazard rates
  hazard_rates <- ifelse(person_time > 0, n_events / person_time, 0)

  # Create hazard table
  hazard_table <- data.frame(
    interval = interval_labels,
    n_at_risk_start = n_risk[1:n_intervals],
    n_censored_interval = n_censored_interval,
    n_events = n_events,
    hazard_rate = hazard_rates,
    stringsAsFactors = FALSE
  )

  # Calculate median survival time
  median_survival <- getMediansurv(time_points, hazard_rates)

  # Return results as list
  return(list(
    hazard_table = hazard_table,
    median_survival = median_survival
  ))
}
