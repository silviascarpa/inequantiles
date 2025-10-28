#' Quintile Share Ratio Estimator
#'
#'
#' Computes the quintile share ratio (QSR) estimator for measuring inequality
#' on simple and complex sampling data
#'
#' @param y A numeric vector of data values
#' @param weights A numeric vector of sampling weights (optional)
#' @param type Quantile estimation type: integer 4-9 or HD for Harrell-Davis (default: 6).
#'         See \code{csquantile} documentation for a complete description.
#' @param na.rm Logical, should missing values be removed? (default: TRUE)
#'
#' @return A scalar numeric value representing the estimated inequality by the
#' quintile share ratio (QSR).
#'
#' @details
#' Consider a random sample \eqn{s} of size \eqn{n}, and let \eqn{w_j}, \eqn{j \in s}, define
#' the sampling weight and \eqn{y_j} be the observed characteristics (i.e. income)
#' associated to the \eqn{j}-th individual, \eqn{j = 1, \ldots, n}.
#' The QSR estimator is defined as:
#'
#' \deqn{\widehat{QSR} =
#' \frac{\sum_{j \in s}w_j y_j \mathds{1}\left\{ y_j \geq \widehat{Q}(0.8)\right\} }{\sum_{j \in s} w_j y_j\mathds{1}\left\{ y_j \leq \widehat{Q}(0.2)\right\} }}
#'
#' where the estimated quantiles \eqn{\widehat{Q}(p)} are computed via the function
#' \code{csquantile()}, which accounts for sampling weights and the specified
#' quantile type. This allows \eqn{\widehat{QSR}} to be used both for simple
#' random samples and for complex survey data with design weights.
#'
#' See  \insertCite{langel2011quintile;textual}{inequantiles} for a complete review of the
#' QSR estimator with complex sampling data.
#'
#'
#'
#' @examples
#'
#' data(synthouse)
#'
#' eq <- synthouse$eq_income ### Income data
#'
#' # Compute unweighted QSR with default type 6 quantile estimator
#' qsr(y = eq)
#'
#' # Consider the sampling weights and change quantile estimation type
#' w <- synthouse$weight
#' qsr(y = eq, weights = w, type = 5)
#'
#' # Compare QSR across macro-regions (NUTS1)
#' tapply(1:nrow(synthouse), synthouse$NUTS1, function(area) {
#'   qsr(y = synthouse$eq_income[area],
#'       weights = synthouse$weight[area],
#'       type = 6)
#' })
#'
#' @references
#'
#' \insertRef{langel2011quintile}{inequantiles}
#'
#' @importFrom Rdpack reprompt
#'
#'
#' @export

### Quintile share ratio
qsr <- function(y, weights = NULL, type = 6, na.rm = TRUE) {

  # Set weights to 1 if not provided
  if (is.null(weights)) {
    weights <- rep(1, length(y))
  }


  # Validate that y and weights have the same length
  if (length(y) != length(weights)) {
    stop("'y' and 'weights' must have the same length")
  }

  # Compute quantiles (csquantile handles na.rm internally)
  q80 <- unname(csquantile(y, 0.8, weights = weights, type = type, na.rm = na.rm))
  q20 <- unname(csquantile(y, 0.2, weights = weights, type = type, na.rm = na.rm))

  # Calculate QSR with weights
  numerator <- sum(weights * y * ifelse(y >= q80, 1, 0), na.rm = na.rm)
  denominator <- sum(weights * y * ifelse(y <= q20, 1, 0), na.rm = na.rm)

  # Check for zero denominator
  if (denominator == 0) {
    warning("Denominator is zero: no observations in the bottom 20%. Returning NA.")
    return(NA_real_)
  }

  est <- numerator / denominator
  return(est)
}

