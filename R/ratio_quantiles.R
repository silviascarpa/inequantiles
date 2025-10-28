#' Ratio of Quantiles Estimator (e.g., P90/P10)
#'
#'
#' Estimates ratio of quantiles (e.g., P90/P10)
#' on simple and complex sampling data
#'
#' @param y A numeric vector of data values
#' @param weights A numeric vector of sampling weights (optional)
#' @param prob_top The percentile to be considered at the numerator (default \code{0.90})
#' @param prob_bottom The percentile to be considered at the denominator (default \code{0.10}
#' @param type Quantile estimation type: integer 4-9 or HD for Harrell-Davis (default: 6).
#'         See \code{csquantile} documentation for a complete description.
#' @param na.rm Logical, should missing values be removed? (default: TRUE)
#'
#' @return A scalar numeric value representing the estimated ratio of quantiles
#'
#' @details
#' Consider a random sample \eqn{s} of size \eqn{n}, and let \eqn{w_j}, \eqn{j \in s}, define
#' the sampling weight and \eqn{y_j} be the observed characteristics (i.e. income)
#' associated to the \eqn{j}-th individual, \eqn{j = 1, \ldots, n}.
#' Let \eqn{p_{top}} be the order of the quantile at the numerator and
#' \eqn{p_{bottom}} be the order of the quantile at the denominator. For example, set \eqn{p_{top} = 0.90} and
#' \eqn{p_{bottom} = 0.10}. Then the popular P90/P10 ratio can be estimated by
#'
#' \deqn{\widehat{{P}90/{P}10} = \frac{\widehat{Q}(p=0.9)}{\widehat{Q}(p=0.1)}  }
#'
#' where the estimated quantiles \eqn{\widehat{Q}(p)} are computed via the function
#' \code{csquantile()}, which accounts for sampling weights and the specified
#' quantile type.
#'
#'
#'
#' @examples
#'
#' data(synthouse)
#'
#' eq <- synthouse$eq_income ### Income data
#'
#' # Compute unweighted P90/P10 with default type 6 quantile estimator
#' ratio_quantiles(y = eq)
#'
#' # Consider the sampling weights, change quantile estimation type and orders of quantiles
#' w <- synthouse$weight
#' ratio_quantiles(y = eq, weights = w, prob_top = 0.6, prob_bottom = 0.1, type = 5)
#'
#' # Compare the P90/P10 across macro-regions (NUTS1)
#' tapply(1:nrow(synthouse), synthouse$NUTS1, function(area) {
#'   ratio_quantiles(y = synthouse$eq_income[area],
#'       weights = synthouse$weight[area])
#' })
#'
#'
#' @importFrom Rdpack reprompt
#'
#'
#' @export

# ------- P90 / P10 ----------
ratio_quantiles <- function(y, weights = NULL, prob_top = 0.90,
                    prob_bottom = 0.10, type = 6, na.rm = FALSE){

  # Set weights to 1 if not provided
  if (is.null(weights)) {
    weights <- rep(1, length(y))
  }

  # Validate that y and weights have the same length
  if (length(y) != length(weights)) {
    stop("'y' and 'weights' must have the same length")
  }

  numerator <- unname(csquantile(y, probs = prob_top, weights = weights, type = type))
  denominator <- unname(csquantile(y, probs = prob_bottom, weights = weights, type = type))

  # Check for zero denominator
  if (denominator == 0) {
    warning("Denominator is zero: Q(", prob_bottom, ") = 0. Returning NA.")
    return(NA_real_)
  }


  ratio <- numerator / denominator
  return(ratio)
}

