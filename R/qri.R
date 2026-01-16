#' Quantile Ratio Index Estimator
#'
#' Computes the quantile ratio index (QRI) estimator for measuring inequality
#'  on simple and complex sampling data
#'
#' @param y A numeric vector of data values
#' @param weights A numeric vector of sampling weights (optional)
#' @param M Integer, number of quantile ratios to average (default: 100)
#' @param type Quantile estimation type: integer 4-9 or HD for Harrell-Davis (default: 6).
#'         See \code{csquantile} documentation for a complete description.
#' @param na.rm Logical, should missing values be removed? (default: TRUE)
#'
#' @return A scalar numeric value representing the estimated inequality by the quantile ratio index (QRI).
#'
#' @details
#' Consider a random sample \eqn{s} of size, where \eqn{w_j}, \eqn{j \in s}, defines
#' the sampling weight associated to the \eqn{j}-th individual.
#' The QRI estimator is defined as:
#'
#' \deqn{\widehat{QRI} = \frac{1}{M}\sum_{m=1}^M\left(1- \frac{\widehat{Q}(p_{m/2})}{\widehat{Q}(1 - p_{m/2})}\right)}
#'
#' where \eqn{p_m = p_m = (m-0.5)/M}.
#' The estimated quantiles \eqn{\widehat{Q}(p)} are computed via the function
#' \code{csquantile()}, which accounts for sampling weights and the specified
#' quantile type. This allows \eqn{\widehat{QRI}} to be used both for simple
#' random samples and for complex survey data with design weights.
#'
#' This index was proposed by \insertCite{prendergast2018simple;textual}{inequantiles},
#' and extended to survey data by \insertCite{scarpa2025inference;textual}{inequantiles}.
#'
#'
#' @examples
#'
#' data(synthouse)
#'
#' eq <- synthouse$eq_income ### Income data
#'
#' # Compute unweighted QRI with default type 6 quantile estimator
#' qri(y = eq)
#'
#' # Consider the sampling weights and change quantile estimation type
#' w <- synthouse$weight
#' qri(y = eq, weights = w, type = 5)
#'
#' # Compare QRI across macro-regions (NUTS1)
#' tapply(1:nrow(synthouse), synthouse$NUTS1, function(area) {
#'   qri(y = synthouse$eq_income[area],
#'       weights = synthouse$weight[area],
#'       type = 6)
#' })
#'
#' @references
#' \insertRef{prendergast2018simple}{inequantiles}
#'
#' \insertRef{scarpa2025inference}{inequantiles}
#'
#' @importFrom Rdpack reprompt
#'
#'
#' @export


qri <- function(y, weights = NULL, M = 100, type = 6, na.rm = TRUE) {

  # Check for negative values
  if (any(y < 0, na.rm = TRUE)) {
    warning("The input vector 'y' contains negative values. ",
            "The QRI is designed for non-negative distributions. ",
            "Results may not be meaningful for data with negative values.",
            call. = FALSE)
  }

  # Set weights to 1 if not provided
  if (is.null(weights)) {
    weights <- rep(1, length(y))
  }

  # Validate that y and weights have the same length
  if (length(y) != length(weights)) {
    stop("'y' and 'weights' must have the same length")
  }

  # Generate probability points
  probs <- ((1:M) - 0.5) / M

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

