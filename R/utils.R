#' Calculate Median Survival Time
#'
#' Calculates the median survival time for a given piecewise exponential model.
#'
#' @param time_points Vector of time points defining the intervals
#' @param hazard_rates Vector of hazard rates for each interval
#'
#' @return The median survival time
#'
#' @details
#' For a piecewise exponential model, the survival function at time t is:
#' S(t) = exp(-sum(hazard_rates[i] * (min(t, time_points[i+1]) - time_points[i])))
#' The median is the time t where S(t) = 0.5
#'
#' @examples
#' time_points <- c(0, 6, 12, 18, 24)
#' hazard_rates <- c(0.05, 0.08, 0.12, 0.15)
#' median_surv <- getMedianSurv(time_points, hazard_rates)
#' print(median_surv)
#'
#' @export
getMedianSurv <- function(time_points, hazard_rates) {
  # Calculate the median survival time

  # Input validation
  if (length(time_points) != length(hazard_rates) + 1) {
    stop("Length of time_points should be length of hazard_rates + 1")
  }

  # Function to calculate survival at time t
  survFun <- function(t) {
    # Accumulate hazard
    total_hazard <- 0

    for (i in seq_along(hazard_rates)) {
      if (t <= time_points[i]) {
        # t is before this interval starts
        break
      }

      # Add hazard for this interval
      interval_length <- min(t, time_points[i+1]) - time_points[i]
      total_hazard <- total_hazard + hazard_rates[i] * interval_length
    }

    return(exp(-total_hazard))
  }

  # Function to find the time where survival = 0.5
  medianFun <- function(t) {
    return(survFun(t) - 0.5)
  }

  # If hazard is 0 everywhere, return Inf
  if (all(hazard_rates == 0)) {
    return(Inf)
  }

  # Otherwise, find the root
  max_time <- time_points[length(time_points)] +
    (time_points[length(time_points)] - time_points[1])

  # Check if survival at max_time is > 0.5
  if (survFun(max_time) > 0.5) {
    # Need to look further
    while (survFun(max_time) > 0.5) {
      max_time <- max_time * 2
    }
  }

  # Find root
  tryCatch({
    median <- stats::uniroot(medianFun, c(0, max_time))$root
    return(median)
  }, error = function(e) {
    warning("Could not find median survival time: ", e$message)
    return(NA)
  })
}

#' Simulate Survival Data from Piecewise Exponential Distribution
#'
#' Generates survival data from a piecewise exponential distribution with
#' specified hazard rates at different time intervals.
#'
#' @param n Number of subjects to simulate
#' @param time_points Vector of time points defining the intervals
#' @param hazard_rates Vector of hazard rates for each interval
#' @param max_time Maximum follow-up time (censoring time)
#' @param censoring_prob Probability of random censoring
#'
#' @return A data frame with columns:
#'   \item{id}{Subject ID}
#'   \item{time}{Event or censoring time}
#'   \item{status}{Event indicator (1=event, 0=censored)}
#'
#' @details
#' The function simulates survival times from a piecewise exponential distribution
#' with specified hazard rates at different time intervals. Random censoring
#' is applied based on the specified probability.
#'
#' @examples
#' time_points <- c(0, 6, 12, 18, 24)
#' hazard_rates <- c(0.05, 0.08, 0.12, 0.15)
#' data <- simPE(
#'   n = 100,
#'   time_points = time_points,
#'   hazard_rates = hazard_rates,
#'   max_time = 36,
#'   censoring_prob = 0.1
#' )
#' head(data)
#'
#' @export
simPE <- function(n, time_points, hazard_rates, max_time = NULL, censoring_prob = 0) {
  # Input validation
  if (length(time_points) != length(hazard_rates) + 1) {
    stop("Length of time_points should be length of hazard_rates + 1")
  }

  if (is.null(max_time)) {
    max_time <- max(time_points)
  }

  # Generate survival times for each subject
  survival_times <- numeric(n)

  for (i in 1:n) {
    # Generate a random number from Uniform(0,1)
    u <- stats::runif(1)

    # Convert to survival time
    t <- 0
    S <- 1

    for (j in seq_along(hazard_rates)) {
      # Calculate survival if event happens in this interval
      lambda <- hazard_rates[j]
      t_start <- time_points[j]
      t_end <- time_points[j+1]

      # Skip if hazard is 0
      if (lambda == 0) {
        t <- t_end
        next
      }

      # Calculate time-to-event within this interval
      interval_surv <- exp(-lambda * (t_end - t_start))

      if (S * interval_surv <= u) {
        # Event occurs in this interval
        t <- t_start - log(S/u) / lambda
        break
      } else {
        # Event occurs in a later interval
        S <- S * interval_surv
        t <- t_end
      }

      # Check if we've reached the end
      if (j == length(hazard_rates)) {
        # Event occurs after last interval
        if (S > u) {
          # Need to extend beyond the last interval
          t <- t_end - log(S/u) / lambda
        }
      }
    }

    survival_times[i] <- t
  }

  # Apply administrative censoring at max_time
  status <- ifelse(survival_times <= max_time, 1, 0)
  times <- pmin(survival_times, max_time)

  # Apply random censoring
  if (censoring_prob > 0) {
    censor_indicator <- stats::rbinom(n, 1, censoring_prob)
    random_censor_times <- stats::runif(n, 0, max_time)

    # Apply random censoring where indicated
    for (i in 1:n) {
      if (censor_indicator[i] == 1 && random_censor_times[i] < times[i]) {
        times[i] <- random_censor_times[i]
        status[i] <- 0
      }
    }
  }

  # Create data frame
  result <- data.frame(
    id = 1:n,
    time = times,
    status = status
  )

  return(result)
}

#' Create Kaplan-Meier Summary Table
#'
#' Creates a summary table of at-risk counts, events, and censored observations
#' at specified time points from survival data.
#'
#' @param data A data frame with columns for time and status
#' @param time_points Vector of time points at which to evaluate
#'
#' @return A list with components:
#'   \item{n_risk}{Number at risk at each time point}
#'   \item{n_event}{Number of events in each interval}
#'   \item{n_censored}{Number censored in each interval}
#'   \item{n_censored_cumulative}{Cumulative number censored up to each time point}
#'
#' @details
#' The function evaluates the number at risk, events, and censored observations
#' at each time point specified in `time_points`. The resulting table is similar
#' to what would be published with Kaplan-Meier curves.
#'
#' @examples
#' time_points <- c(0, 6, 12, 18, 24)
#' hazard_rates <- c(0.05, 0.08, 0.12, 0.15)
#' data <- simPE(
#'   n = 100,
#'   time_points = time_points,
#'   hazard_rates = hazard_rates,
#'   max_time = 36,
#'   censoring_prob = 0.1
#' )
#' km_table <- summaryKM(data, time_points)
#' print(km_table)
#'
#' @export
summaryKM <- function(data, time_points) {
  # Ensure time_points is sorted
  time_points <- sort(time_points)

  # Extract data
  times <- data$time
  status <- data$status

  # Initialize vectors
  n_risk <- numeric(length(time_points))
  n_event <- numeric(length(time_points) - 1)
  n_censored <- numeric(length(time_points) - 1)
  n_censored_cumulative <- numeric(length(time_points))

  # Calculate n_risk at each time point
  for (i in seq_along(time_points)) {
    n_risk[i] <- sum(times >= time_points[i])
  }

  # Calculate events and censoring in each interval
  for (i in 1:(length(time_points) - 1)) {
    t_start <- time_points[i]
    t_end <- time_points[i+1]

    # Events in this interval
    n_event[i] <- sum(times >= t_start & times < t_end & status == 1)

    # Censored in this interval
    n_censored[i] <- sum(times >= t_start & times < t_end & status == 0)
  }

  # Calculate cumulative censoring
  for (i in 1:length(time_points)) {
    if (i == 1) {
      n_censored_cumulative[i] <- 0
    } else {
      n_censored_cumulative[i] <- sum(n_censored[1:(i-1)])
    }
  }

  # Return results as a list
  return(list(
    n_risk = n_risk,
    n_event = n_event,
    n_censored = n_censored,
    n_censored_cumulative = n_censored_cumulative
  ))
}

#' Extract Hazard Rates from Kaplan-Meier Data
#'
#' Estimates hazard rates from Kaplan-Meier data using the person-years approach.
#'
#' @param time_points Vector of time points defining the intervals
#' @param n_risk Vector of number at risk at each time point
#' @param n_censored Vector of cumulative number censored up to each time point
#' @param n_events Optional vector of events in each interval; if NULL, calculated from n_risk
#' @param warn_negative_events Logical, whether to warn about negative event counts
#'
#' @return A list with components:
#'   \item{hazard_table}{Data frame with estimated hazard rates}
#'   \item{median_survival}{Estimated median survival time}
#'   \item{survival_function}{Function to calculate survival at any time}
#'
#' @details
#' The function estimates hazard rates from Kaplan-Meier data using the person-years
#' approach. It calculates the number of events in each interval (if not provided),
#' the person-time at risk, and the hazard rate.
#'
#' @examples
#' time_points <- c(0, 6, 12, 18, 24)
#' hazard_rates <- c(0.05, 0.08, 0.12, 0.15)
#' data <- simPE(
#'   n = 100,
#'   time_points = time_points,
#
