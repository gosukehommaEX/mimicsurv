#' Create Kaplan-Meier Style Summary Table from Simulated Data
#'
#' This function creates a summary table in Kaplan-Meier format from simulated
#' individual patient data, calculating number at risk and cumulative censoring
#' at specified time points.
#'
#' @param sim_data Data frame with time and status columns from simulation
#' @param time_points Time points for summary table creation
#'
#' @return List containing:
#'   \item{n_risk}{Numeric vector of number at risk at each time point}
#'   \item{n_censored_cumulative}{Numeric vector of cumulative censored by each time point}
#'
#' @details
#' The function processes individual patient data to create aggregate summaries
#' that mimic published Kaplan-Meier tables. Number at risk represents patients
#' who have not yet experienced the event or been censored before each time point.
#'
#' @examples
#' # First simulate data
#' time_points <- c(0, 6, 12, 18, 24)
#' hazard_rates <- c(0.05, 0.08, 0.12, 0.15)
#'
#' sim_data <- simPE(
#'   n = 200,
#'   time_points = time_points,
#'   hazard_rates = hazard_rates,
#'   max_time = 24
#' )
#'
#' # Create KM table
#' km_summary <- summaryKM(sim_data, time_points)
#' print(km_summary$n_risk)
#' print(km_summary$n_censored_cumulative)
#'
#' @export

summaryKM <- function(sim_data, time_points) {

  # Input validation
  if (!is.data.frame(sim_data)) {
    stop("sim_data must be a data frame")
  }

  if (!all(c("time", "status") %in% names(sim_data))) {
    stop("sim_data must contain 'time' and 'status' columns")
  }

  n <- nrow(sim_data)
  n_timepoints <- length(time_points)

  n_risk <- numeric(n_timepoints)
  n_censored_cumulative <- numeric(n_timepoints)

  for (i in 1:n_timepoints) {
    t <- time_points[i]

    # Number at risk: not yet had event or censored before time t
    n_risk[i] <- sum(sim_data$time >= t)

    # Cumulative number censored by time t
    n_censored_cumulative[i] <- sum(sim_data$time <= t & sim_data$status == 0)
  }

  return(list(
    n_risk = n_risk,
    n_censored_cumulative = n_censored_cumulative
  ))
}
