#' Ratio of Quantiles (e.g., P90/P10)
#'
#'
#' Estimates ratio of quantiles (e.g., P90/P10) on simple and complex sampling data
#'
#' @param y A numeric vector of data values
#' @param weights A numeric vector of sampling weights (optional). If \code{NULL},
#'   all observations are equally weighted.
#' @param prob_numerator The percentile to be considered at the numerator (default \code{0.90})
#' @param prob_denominator The percentile to be considered at the denominator (default \code{0.10})
#' @param type Quantile estimation type: integer \code{4}--\code{9} or
#'   \code{"HD"} for Harrell--Davis (default: \code{6}). See \code{\link{csquantile}}.
#' @param na.rm Logical, should missing values be removed? (default: TRUE)
#'
#' @returns A scalar numeric value representing the estimated ratio of quantiles
#'
#' @details
#' Consider a random sample \eqn{s} of size \eqn{n}, and let \eqn{w_j}, \eqn{j \in s}, define
#' the sampling weight and \eqn{y_j} be the observed characteristics (i.e. income)
#' associated to the \eqn{j}-th individual, \eqn{j = 1, \ldots, n}.
#' Let \eqn{p_{n}} be the order of the quantile at the numerator and
#' \eqn{p_{d}} be the order of the quantile at the denominator. For example, set \eqn{p_{n} = 0.90} and
#' \eqn{p_{d} = 0.10}. Then the popular P90/P10 ratio can be estimated by
#'
#' \deqn{\widehat{{P}90/{P}10} = \frac{\widehat{Q}(p_n=0.9)}{\widehat{Q}(p_d=0.1)}  }
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
#' eq <- synthouse$eq_income # Income data
#'
#' # Compute unweighted P90/P10 with default type 6 quantile estimator
#' ratio_quantiles(y = eq)
#'
#' # Consider the sampling weights, change quantile estimation type and orders of quantiles
#' w <- synthouse$weight
#' ratio_quantiles(y = eq, weights = w, prob_numerator = 0.6, prob_denominator = 0.1, type = 5)
#'
#' # Compare the P90/P10 across macro-regions (NUTS1)
#' tapply(1:nrow(synthouse), synthouse$NUTS1, function(area) {
#'   ratio_quantiles(y = synthouse$eq_income[area],
#'       weights = synthouse$weight[area])
#' })
#'
#'
#' @family inequality indicators based on quantiles
#'
#' @importFrom Rdpack reprompt
#'
#'
#' @export

# ------- P90 / P10 ----------
ratio_quantiles <- function(y, weights = NULL, prob_numerator = 0.90,
                    prob_denominator = 0.10, type = 6, na.rm = TRUE){

  # Set weights to 1 if not provided
  if (is.null(weights)) {
    weights <- rep(1, length(y))
  }

  # Validate that y and weights have the same length
  if (length(y) != length(weights)) {
    stop("'y' and 'weights' must have the same length")
  }

  numerator <- unname(csquantile(y, probs = prob_numerator, weights = weights,
                                 type = type, na.rm = na.rm))
  denominator <- unname(csquantile(y, probs = prob_denominator, weights = weights,
                                   type = type, na.rm = na.rm))
  # Check for zero denominator
  if (denominator == 0) {
    warning("Denominator is zero: Q(", prob_denominator, ") = 0. Returning NA.")
    return(NA_real_)
  }


  ratio <- numerator / denominator

  name_num <- prob_numerator * 100
  name_denom <- prob_denominator * 100
  names(ratio) <- paste0("P", name_num, "/", "P", name_denom)

  return(ratio)
}

