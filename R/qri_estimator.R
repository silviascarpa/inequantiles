#' quantile ratio index estimator
#'
#' Computes the quantile ratio index (QRI) estimator for measuring inequality
#'
#' @param y A numeric vector of data values
#' @param weights A numeric vector of sampling weights (optional)
#' @param J Integer, number of quantile ratios to average (default: 100)
#' @param type Quantile estimation type: integer 4-9 or "HD" for Harrell-Davis (default: 6)
#' @param na.rm Logical, should missing values be removed? (default: TRUE)
#'
#' @return A scalar numeric value representing the estimated quantile ratio index (QRI)
#'
#' @details
#' The QRI estimator is defined as:
#'
#' \deqn{\widehat{QRI} = \frac{1}{M}\sum_{m=1}^M\left(1- \frac{\widehat{Q}(p_{m/2})}{\widehat{Q}(1 - p_{m/2})}\right)}
#'
#' where the estimated quantiles \eqn{\widehat{Q}(p)} are computed via the function
#' \code{csquantile()}, which accounts for sampling weights and the specified
#' quantile type. This allows \eqn{\widehat{QRI}} to be used both for simple
#' random samples and for complex survey data with design weights.
#'
#'
#' @references
#' Scarpa, S., Ferrante, M. R., & Sperlich, S. (2025).
#' "Inference for the Quantile Ratio Inequality Index in the Context of Survey Data."
#' *Journal of Survey Statistics and Methodology*.
#'
#'
#' @export


qri <- function(y, weights = NULL, J = 100, type = 6, na.rm = TRUE) {

  # Generate probability points
  probs <- ((1:J) - 0.5) / J

  # Compute quantiles at p/2 and 1-p/2
  q_lower <- csquantile(y, weights = weights, probs = probs / 2,
                        type = type, na.rm = na.rm)
  q_upper <- csquantile(y, weights = weights, probs = 1 - probs / 2,
                        type = type, na.rm = na.rm)

  # Compute ratios
  Rp <- q_lower / q_upper

  # Compute QRI
  qri_value <- mean(1 - Rp)

  return(qri_value)
}

