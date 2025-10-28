#' Quantile Estimator for Grouped (Binned) Data
#'
#' Computes quantiles from grouped frequency data using linear interpolation
#' within the quantile class.
#'
#' @param freq Numeric vector of class frequencies (counts). Must be non-negative.
#' @param lower_bounds Numeric vector of lower class bounds. Must be strictly
#'   increasing.
#' @param upper_bounds Numeric vector of upper class bounds. Must be strictly
#'   increasing and greater than corresponding \code{lower_bounds}.
#' @param probs Numeric vector of probabilities (between 0 and 1) for which to
#'   compute the quantiles. Default is 0.5 (median).
#' @param midpoints Optional numeric vector of class midpoints. Used only as
#'   fallback when a quantile class has zero frequency. If \code{NULL}, the
#'   midpoint is computed as the arithmetic mean of class bounds.
#'
#' @return A vector of estimated quantiles on grouped data corresponding to \code{probs}.
#'   Returns \code{NA} if total frequency is zero or missing.
#'
#' @details
#'
#' Consider grouped data divided into \eqn{L} classes with known boundaries. Let:
#' \itemize{
#'   \item \eqn{L_j} be the lower bound of the \eqn{j}-th quantile class
#'   \item \eqn{U_j} be the upper bound of the \eqn{j}-th quantile class
#'   \item \eqn{h_j = U_j - L_j} be the \eqn{j}-th quantile class width
#'   \item \eqn{C_{j-1}} be the cumulative frequency up to the previous class
#'   \item \eqn{f_j} be the frequency within the quantile class \eqn{j}
#'   \item \eqn{N = \sum_{i=1}^{k} f_i} be the total frequency
#' }
#'
#' The quantile class for the \eqn{p}-th quantile is the first class \eqn{j} such that:
#'
#' \deqn{j = min\{i: C_i \geq pN \}}.
#'
#' The \eqn{p}-th quantile \eqn{Q(p)} is then estimated by linear interpolation within the
#' quantile class:
#'
#' \deqn{\widetilde{Q(p)} = L_j + \frac{(pN - C_{j-1})}{f_j} \cdot h_j}
#'
#'
#' The method assumes a uniform distribution of observations within
#' each class interval. This is a standard approach for grouped data when individual
#' observations are not available.
#'
#'@section Handling Open-Ended Classes:
#'
#' When dealing with administrative or tax data, the first class often is often defined
#' as negative income  (or incomes below zero) and the last class as incomes above
#' a certain threshold. In such cases, we have \code{-Inf} as the lower bound of the
#' first class and \code{Inf} as the upper bound of the last class.
#'
#' If \code{Inf} values are present in the given bounds, the function imputes reasonable
#' bounds using the specified method:
#'
#' \strong{For open left class (first lower bound = -Inf):}
#' The imputed first lower bound is given by:
#' \deqn{L_1^* = U_1 - h_2}
#'
#' where \eqn{U_1} is the upper bound of the first class and \eqn{h_2 = U_2 - L_2}
#' is the width of the second class. This assumes the first class has the same
#' width as the second class.
#'
#' \strong{For open right class (last upper bound = Inf):}
#'
#' The imputed upper bound is given by:
#' \deqn{U_J^* = L_J + h_{J-1}}
#'
#' where \eqn{L_J} is the lower bound of the last class and
#' \eqn{h_{J-1} = U_{J-1} - L_{J-1}} is the width of the second-to-last class.
#' This assumes the last class has the same width as the penultimate class.
#'
#' @section Special Cases:
#' \itemize{
#'   \item If the quantile class has zero frequency, the function returns the
#'     class midpoint as a fallback.
#'   \item If total frequency is zero or \code{NA}, the function returns \code{NA}
#'     for all requested quantiles.
#' }
#'
#'
#' @examples
#' # Basic usage: compute quartiles
#' freq <- c(5, 8, 10, 4, 3)
#' lower <- c(0, 10, 20, 30, 40)
#' upper <- c(10, 20, 30, 40, 50)
#'
#' quantile_grouped(freq, lower, upper, probs = c(0.25, 0.5, 0.75))
#'
#' # Compute deciles
#' quantile_grouped(freq, lower, upper, probs = seq(0.1, 0.9, by = 0.1))
#'
#' # With custom midpoints
#' midpts <- c(5, 15, 25, 35, 45)
#' quantile_grouped(freq, lower, upper, probs = 0.5, midpoints = midpts)
#'
#' # Income distribution example
#' income_freq <- c(120, 180, 150, 80, 40, 20, 10)
#' income_lower <- c(0, 15000, 30000, 45000, 60000, 80000, 100000)
#' income_upper <- c(15000, 30000, 45000, 60000, 80000, 100000, 150000)
#'
#' # Compute median income
#' quantile_grouped(income_freq, income_lower, income_upper, probs = 0.5)
#'
#' # Compute income quintiles
#' quantile_grouped(income_freq, income_lower, income_upper,
#'                    probs = seq(0.2, 0.8, by = 0.2))
#'
#'
#'
#' @seealso
#' \code{\link{quantile}} for quantiles of ungrouped data.
#'
#' @family grouped data functions
#'
#' @export
quantile_grouped <- function(freq, lower_bounds, upper_bounds,
                               probs = 0.5, midpoints = NULL) {

  # Input validation
  if (length(freq) != length(lower_bounds) || length(freq) != length(upper_bounds)) {
    stop("freq, lower_bounds, and upper_bounds must have the same length")
  }

  if (any(probs < 0 | probs > 1, na.rm = TRUE)) {
    stop("probs must be between 0 and 1")
  }

  if (any(freq < 0, na.rm = TRUE)) {
    warning("Negative frequencies detected and will be treated as NA")
  }

  if (!is.null(midpoints) && length(midpoints) != length(freq)) {
    stop("midpoints must have the same length as freq")
  }

  # Calculate total frequency
  freq <- as.numeric(freq)
  total <- sum(freq, na.rm = TRUE)

  # Return NA for all quantiles if total is zero or NA
  if (is.na(total) || total == 0) {
    return(rep(NA_real_, length(probs)))
  }

  # Handle infinite bounds (replace with finite approximations)
  if (is.infinite(lower_bounds[1])) {
    lower_bounds[1] <- upper_bounds[1] - (upper_bounds[2] - lower_bounds[2])
  }
  if (is.infinite(upper_bounds[length(upper_bounds)])) {
    upper_bounds[length(upper_bounds)] <- lower_bounds[length(lower_bounds)] +
      (lower_bounds[length(lower_bounds)] - lower_bounds[length(lower_bounds)-1])
  }

  # Compute cumulative frequencies
  cum <- cumsum(freq)

  # Initialize result vector
  res <- numeric(length(probs))

  # Compute each quantile
  for (i in seq_along(probs)) {
    p <- probs[i]
    pos <- total * p

    # Find the quantile class
    idx <- which(cum >= pos)[1]

    # Fallback to last class if not found
    if (is.na(idx)) {
      idx <- length(freq)
    }

    # Extract class parameters
    L <- lower_bounds[idx]
    CumPrev <- if (idx == 1) 0 else cum[idx - 1]
    FreqClass <- freq[idx]
    h <- upper_bounds[idx] - lower_bounds[idx]

    # Handle zero-frequency class
    if (is.na(FreqClass) || FreqClass == 0) {
      if (!is.null(midpoints)) {
        res[i] <- midpoints[idx]
      } else {
        mid <- (lower_bounds + upper_bounds) / 2
        res[i] <- mid
      }
    } else {
      # Linear interpolation within the class
      res[i] <- L + ((pos - CumPrev) / FreqClass) * h
    }
  }

  res <- as.numeric(res)
  names(res) <- paste(probs * 100, "%")

  return(res)
}


