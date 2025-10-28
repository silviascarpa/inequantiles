#' Quantile Ratio Rndex Estimator for Grouped (Binned) Data
#'
#' Computes the quantile ratio index (QRI) for measuring inequality from
#' grouped frequency data using linear interpolation for quantile estimation.
#' This function is intended to be used for administrative or tax data, which are
#' very often in the form of grouped (binned) data. Therefore, sampling weights are
#' not considered.
#'
#' @param freq Numeric vector of class frequencies (counts). Must be non-negative.
#' @param lower_bounds Numeric vector of lower class bounds.
#' @param upper_bounds Numeric vector of upper class bounds.
#' @param J Integer, number of quantile ratios to average (default: 100).
#' @param midpoints Optional numeric vector of class midpoints. Used only as
#'   fallback when a quantile class has zero frequency.
#' @param na.rm Logical, should missing values in frequencies be removed? (default: TRUE)
#'
#' @return A scalar numeric value representing the estimated inequality by the
#'   quantile ratio index (QRI) for grouped data.
#'
#' @details
#' Consider grouped data divided into \eqn{L} classes with known boundaries and
#' observed frequencies \eqn{f_1, \ldots, f_L}. The QRI estimator for grouped
#' data is defined as:
#'
#' \deqn{\widetilde{QRI}_{\text{grouped}} = \frac{1}{J}\sum_{j=1}^{J}\left(1 - \frac{\widetilde{Q}(p_j/2)}{\widetilde{Q}(1 - p_j/2)}\right)}
#'
#' where:
#' \itemize{
#'   \item \eqn{p_j = (j - 0.5)/J} for \eqn{j = 1, \ldots, J}
#'   \item \eqn{\widetilde{Q}(p)} denotes the \eqn{p}-th quantile computed from
#'     grouped data using linear interpolation (see \code{\link{quantile_grouped}})
#'   \item \eqn{J} is the number of quantile ratios to average (default: 100)
#' }
#'
#' The quantiles \eqn{\widetilde{Q}(p)} are computed via \code{quantile_grouped()},
#' which uses linear interpolation within classes and automatically handles
#' open-ended classes (with \code{-Inf} or \code{Inf} bounds).
#'
#' The QRI ranges from 0 (perfect equality) to 1 (maximum inequality). The index
#' measures inequality by averaging the relative differences between symmetric quantiles
#' below and above the median, across the entire distribution.
#'
#'
#'
#' @section Comparison with Microdata QRI:
#'
#' When individual-level (microdata) are available, use \code{\link{qri}} instead,
#' which provides more accurate estimates. The grouped data version
#' \code{qri_grouped} should be used when only frequency distributions are available,
#' such as in published statistical tables or administrative aggregates.
#'
#' The grouped QRI will generally approximate the microdata QRI well when:
#' \itemize{
#'   \item Classes are sufficiently narrow
#'   \item The distribution within classes is approximately uniform
#'   \item Sample sizes within classes are adequate
#' }
#'
#' @examples
#' # Basic example with closed classes
#' income_freq <- c(120, 180, 150, 80, 40, 20, 10)
#' income_lower <- c(0, 15000, 30000, 45000, 60000, 80000, 100000)
#' income_upper <- c(15000, 30000, 45000, 60000, 80000, 100000, 150000)
#'
#' qri_grouped(income_freq, income_lower, income_upper)
#'
#' # Example with open-ended classes (Italian MEF-style data)
#' wage_freq <- c(150, 200, 180, 220, 180, 50, 15, 5)
#' wage_lower <- c(-Inf, 0, 10000, 15000, 26000, 55000, 75000, 120000)
#' wage_upper <- c(0, 10000, 15000, 26000, 55000, 75000, 120000, Inf)
#'
#' # Compute QRI (automatically handles open classes)
#' qri_grouped(wage_freq, wage_lower, wage_upper)
#'
#'
#'
#' @seealso
#' \code{\link{qri}} for QRI estimation with microdata.
#' \code{\link{quantile_grouped}} for quantile estimation from grouped data.
#'
#' @family inequality measures
#' @family grouped data functions
#'
#' @references
#' \insertRef{prendergast2018simple}{inequantiles}
#'
#' @export
qri_grouped <- function(freq, lower_bounds, upper_bounds,
                        J = 100, midpoints = NULL, na.rm = TRUE) {

  # Input validation
  if (length(freq) != length(lower_bounds) || length(freq) != length(upper_bounds)) {
    stop("freq, lower_bounds, and upper_bounds must have the same length")
  }

  if (!is.numeric(J) || length(J) != 1 || J < 1) {
    stop("J must be a positive integer")
  }

  if (any(freq < 0, na.rm = TRUE)) {
    warning("Negative frequencies detected and will be treated as NA")
  }

  # Check for negative values in bounds (would indicate negative incomes/wealth)
  if (any(lower_bounds[!is.infinite(lower_bounds)] < 0, na.rm = TRUE) ||
      any(upper_bounds[!is.infinite(upper_bounds)] < 0, na.rm = TRUE)) {
    warning("Class bounds contain negative values. ",
            "The QRI is designed for non-negative distributions. ",
            "Results may not be meaningful for data with negative values.",
            call. = FALSE)
  }

  # Handle NA removal
  if (na.rm) {
    valid_idx <- !is.na(freq) & !is.na(lower_bounds) & !is.na(upper_bounds)
    freq <- freq[valid_idx]
    lower_bounds <- lower_bounds[valid_idx]
    upper_bounds <- upper_bounds[valid_idx]
  }

  # Check if there's any data left
  total_freq <- sum(freq, na.rm = TRUE)
  if (is.na(total_freq) || total_freq == 0) {
    warning("Total frequency is zero or NA. Returning NA.")
    return(NA_real_)
  }

  # Generate probability points: p_j = (m - 0.5) / J for m = 1, ..., J
  probs <- ((1:J) - 0.5) / J

  # Compute lower quantiles: Q(p_j / 2)
  q_lower <- quantile_grouped(
    freq = freq,
    lower_bounds = lower_bounds,
    upper_bounds = upper_bounds,
    probs = probs / 2,
    midpoints = midpoints
  )

  # Compute upper quantiles: Q(1 - p_j / 2)
  q_upper <- quantile_grouped(
    freq = freq,
    lower_bounds = lower_bounds,
    upper_bounds = upper_bounds,
    probs = 1 - probs / 2,
    midpoints = midpoints
  )

  # Compute quantile ratios: R_j = Q(p_j / 2) / Q(1 - p_j / 2)
  Rp <- q_lower / q_upper

  # Handle potential division issues
  if (any(is.na(Rp)) || any(is.infinite(Rp))) {
    warning("Some quantile ratios are NA or Inf. This may indicate issues with ",
            "the class boundaries or frequencies.")
  }

  # Compute QRI: average of (1 - R_j)
  qri_value <- mean(1 - Rp, na.rm = TRUE)

  return(qri_value)
}


