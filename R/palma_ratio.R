#' Palma Ratio Estimator
#'
#'
#' Computes the Palma ratio estimator for measuring inequality
#' on simple and complex sampling data
#'
#' @param y A numeric vector of data values
#' @param weights A numeric vector of sampling weights (optional)
#' @param type Quantile estimation type: integer 4-9 or HD for Harrell-Davis (default: 6).
#'         See \code{csquantile} documentation for a complete description.
#' @param na.rm Logical, should missing values be removed? (default: TRUE)
#'
#' @return A scalar numeric value representing the estimated inequality by the
#' Palma ratio.
#'
#' @details
#' Consider a random sample \eqn{s} of size \eqn{n}, and let \eqn{w_j}, \eqn{j \in s}, define
#' the sampling weight and \eqn{y_j} be the observed characteristics (i.e. income)
#' associated to the \eqn{j}-th individual, \eqn{j = 1, \ldots, n}.
#' According to \insertCite{cobham2013all;textual}{inequantiles} definition,
#' the Palma index divides the share earned by the richest 10\% by the share of the
#' poorest 40\%. Its estimator is defined as:
#'
#' \deqn{\widehat{Palma} =
#' \frac{\sum_{j \in s}w_j y_j \mathds{1}\left\{ y_j \geq \widehat{Q}(0.9)\right\} }{\sum_{j \in s} w_j y_j\mathds{1}\left\{ y_j \leq \widehat{Q}(0.4)\right\} }}
#'
#' where the estimated quantiles \eqn{\widehat{Q}(p)} are computed via the function
#' \code{csquantile()}, which accounts for sampling weights and the specified
#' quantile type. This allows \eqn{\widehat{Palma}} to be used both for simple
#' random samples and for complex survey data with design weights.
#'
#' See  \insertCite{palma2006globalizing;textual}{inequantiles} and \insertCite{palma2011homogeneous;textual}{inequantiles}
#' for an introduction to the Palma ratio.
#'
#'
#' @examples
#'
#' data(synthouse)
#'
#' eq <- synthouse$eq_income ### Income data
#'
#' # Compute unweighted Palma index with default type 6 quantile estimator
#' palma_ratio(y = eq)
#'
#' # Consider the sampling weights and change quantile estimation type
#' w <- synthouse$weight
#' palma_ratio(y = eq, weights = w, type = 5)
#'
#' # Compare the Palma index across macro-regions (NUTS1)
#' tapply(1:nrow(synthouse), synthouse$NUTS1, function(area) {
#'   palma_ratio(y = synthouse$eq_income[area],
#'       weights = synthouse$weight[area],
#'       type = 6)
#' })
#'
#' @references
#'
#' \insertRef{palma2006globalizing}{inequantiles}
#'
#' \insertRef{palma2011homogeneous}{inequantiles}
#'
#' \insertRef{cobham2013all}{inequantiles}
#'
#' @importFrom Rdpack reprompt
#'
#'
#' @export

### Palma ratio
palma_ratio <- function(y, weights = NULL, type = 6, na.rm = TRUE) {

  # Set weights to 1 if not provided
  if (is.null(weights)) {
    weights <- rep(1, length(y))
  }

  # Validate that y and weights have the same length
  if (length(y) != length(weights)) {
    stop("'y' and 'weights' must have the same length")
  }

  # Compute quantiles (csquantile handles na.rm internally)
  q90 <- unname(csquantile(y, 0.9, weights = weights, type = type, na.rm = na.rm))
  q40 <- unname(csquantile(y, 0.4, weights = weights, type = type, na.rm = na.rm))

  # Calculate QSR with weights
  numerator <- sum(weights * y * ifelse(y >= q90, 1, 0), na.rm = na.rm)
  denominator <- sum(weights * y * ifelse(y <= q40, 1, 0), na.rm = na.rm)

  # Check for zero denominator
  if (denominator == 0) {
    warning("Denominator is zero: no observations in the bottom 40%. Returning NA.")
    return(NA_real_)
  }

  est <- numerator / denominator
  return(est)
}

