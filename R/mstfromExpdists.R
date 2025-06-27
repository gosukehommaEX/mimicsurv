#' Calculate Overall Median Survival Time from Mixed Exponential Distributions
#'
#' This function calculates the overall median survival time for a mixed population
#' composed of multiple subgroups, each following exponential distributions with
#' different median survival times. The overall population follows a mixture of
#' exponential distributions with mixing probabilities based on expected event counts.
#'
#' @param expected_events Numeric vector of expected event counts for each subgroup
#' @param MST_subgroups Numeric vector of median survival times for each subgroup (same length as expected_events)
#' @param search_range Numeric vector of length 2 specifying the search interval for uniroot (default: c(0, 100))
#'
#' @return List containing:
#'   \item{subgroup_summary}{Data frame with subgroup details (subgroup, expected_events, MST)}
#'   \item{MST_overall}{Overall median survival time for the mixed population}
#'
#' @details
#' The function assumes that each subgroup follows an exponential distribution with
#' the specified median survival time. The overall population is modeled as a mixture
#' of exponential distributions with mixing probabilities proportional to expected event counts.
#'
#' Mathematical formulation:
#' \deqn{S(t) = \sum_{j=1}^{k} p_j \exp(-\lambda_j t)}
#' where \eqn{p_j = E_j / \sum E_j} is the mixing probability for subgroup j based on expected events,
#' \eqn{\lambda_j = \ln(2) / MST_j} is the hazard rate for subgroup j, and
#' \eqn{MST_j} is the median survival time for subgroup j.
#'
#' The overall median survival time is found by solving:
#' \deqn{S(t_{median}) = 0.5}
#'
#' @examples
#' # Example 1: Two subgroups with different expected event counts
#' expected_events <- c(75, 120)  # Expected events, not sample sizes
#' MST_subgroups <- c(6, 7.5)
#' result <- mstfromExpdists(expected_events, MST_subgroups)
#' print(result$subgroup_summary)
#' print(paste("Overall MST:", round(result$MST_overall, 2), "months"))
#'
#' # Example 2: Three subgroups
#' expected_events <- c(150, 180, 60)  # Based on expected event rates
#' MST_subgroups <- c(8, 12, 15)
#' result <- mstfromExpdists(expected_events, MST_subgroups)
#' print(result$subgroup_summary)
#' print(paste("Overall MST:", round(result$MST_overall, 2), "months"))
#'
#' @importFrom stats uniroot
#' @export

mstfromExpdists <- function(expected_events, MST_subgroups, search_range = c(0, 100)) {

  # Input validation
  if (length(expected_events) != length(MST_subgroups)) {
    stop("expected_events and MST_subgroups must have the same length")
  }

  if (length(expected_events) < 2) {
    stop("At least 2 subgroups are required")
  }

  if (any(expected_events <= 0)) {
    stop("All expected event counts must be positive")
  }

  if (any(MST_subgroups <= 0)) {
    stop("All median survival times must be positive")
  }

  if (length(search_range) != 2 || search_range[1] >= search_range[2]) {
    stop("search_range must be a vector of length 2 with search_range[1] < search_range[2]")
  }

  # Calculate hazard rates from median survival times
  lambda <- log(2) / MST_subgroups

  # Calculate mixing probabilities based on expected event counts
  total_expected_events <- sum(expected_events)
  mixing_probs <- expected_events / total_expected_events

  # Define survival equation: S(t) - 0.5 = 0
  survival_equation <- function(t) {
    survival_prob <- sum(mixing_probs * exp(-lambda * t))
    return(survival_prob - 0.5)
  }

  # Check if the equation can be solved in the given range
  f_lower <- survival_equation(search_range[1])
  f_upper <- survival_equation(search_range[2])

  if (f_lower * f_upper > 0) {
    if (f_lower > 0) {
      warning("Median survival time may be beyond the upper search range. Consider increasing search_range[2].")
      return(list(
        subgroup_summary = data.frame(
          subgroup = 1:length(expected_events),
          expected_events = expected_events,
          MST = MST_subgroups
        ),
        MST_overall = as.numeric(NA)
      ))
    } else {
      warning("Median survival time may be below the lower search range. Consider decreasing search_range[1].")
      return(list(
        subgroup_summary = data.frame(
          subgroup = 1:length(expected_events),
          expected_events = expected_events,
          MST = MST_subgroups
        ),
        MST_overall = as.numeric(NA)
      ))
    }
  }

  # Solve for median survival time
  tryCatch({
    result_uniroot <- uniroot(survival_equation, interval = search_range)
    MST_overall <- result_uniroot$root
  }, error = function(e) {
    stop("Failed to find median survival time. Error: ", e$message)
  })

  # Create subgroup summary data frame
  subgroup_summary <- data.frame(
    subgroup = 1:length(expected_events),
    expected_events = expected_events,
    MST = MST_subgroups
  )

  # Return results
  return(list(
    subgroup_summary = subgroup_summary,
    MST_overall = MST_overall
  ))
}
