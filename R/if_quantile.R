#' Influence Function for Quantiles
#'
#' Computes the influence function of sample quantiles, allowing for both
#' simple random sampling and complex survey designs with sampling weights, in
#' the context of finite population. See \insertCite{hampel1986robust;textual}{inequantiles}
#' for an explanation of influence function and \insertCite{deville1999variance;textual}{inequantiles}
#' for its definition in finite population theory.
#'
#'
#' @param y A numeric vector of data values
#' @param weights A numeric vector of sampling weights (optional)
#' @param probs A numeric value specifying the probability for the quantile (e.g., 0.5 for median)
#' @param type Quantile estimation type: integer 4-9 or "HD" for Harrell-Davis (default: 6)
#' @param na.rm Logical, should missing values be removed? (default: TRUE)
#'
#' @return A numeric vector containing the estimtaed quantile influence function values for each observation.
#'
#' @details
#' From the definiton in \insertCite{van2000asymptotic;textual}{inequantiles} and \insertCite{osier2009variance}{inequantiles},
#' the population influence function of the quantile \eqn{Q(p)} is defined as:
#'
#' \deqn{IF(Q(p))_i = \frac{p - \mathbf{1}(y_i \leq Q(p))}{f(Q(p)) \, N},}
#'
#' where \eqn{f(Q(p))} is the population density function evaluated at the quantile and
#'    \eqn{N} is the population size.
#'
#' In the sample, this is estimated as:
#'
#' \deqn{\widehat{IF}(Q(p))_i = \frac{p - \mathbf{1}(y_i \leq \widehat{Q}(p))}{\widehat{f}(\widehat{Q}(p)) \, \widehat{N}},}
#'
#' where \eqn{\widehat{Q}(p)} is the weighted sample quantile estimated by
#' \code{csquantile()}, and \eqn{\widehat{N} = \sum_{i \in s} w_i} is the estimated population size.
#'
#' The density \eqn{\widehat{f}(y)} is estimated using a Gaussian kernel density function:
#'
#' \deqn{
#' \widehat{f}(y) = \frac{1}{\widehat{N} \, h \sqrt{2\pi}}
#' \sum_{j \in s} w_j \exp\!\left\{ -\frac{(y - y_j)^2}{2h^2} \right\},
#' }
#'
#' with bandwidth \eqn{h = 0.79 \cdot IQR \cdot \widehat{N}^{-1/5}}
#'
#'
#' @examples
#'
#' # On synthetic data
#' eq_synth <- rlnorm(30, 9, 0.7)
#' IF_synth <- if_quantile(y = eq_synth, probs = 0.3)
#'
#' # On real data
#' data(synthouse)
#' eq <- synthouse$eq_income[1:30] ## Take some observations (as example)
#' w <- synthouse$weight[1:30]
#' IF_quantile <- if_quantile(y = eq, weights = w, type = 6, probs = 0.5)
#'
#'
#'
#' @references
#' \insertRef{hampel1986robust}{inequantiles}
#'
#' \insertRef{deville1999variance}{inequantiles}
#'
#' \insertRef{van2000asymptotic}{inequantiles}
#'
#' \insertRef{osier2009variance}{inequantiles}
#
#
#'
#'
#' @export
if_quantile <- function(y, weights = NULL, probs, type = 6, na.rm = TRUE) {

  # Compute estimated population size
  if (is.null(weights)) {
    hatN <- length(y)
    weights_internal <- rep(1, length(y))
  } else {
    hatN <- sum(weights)
    weights_internal <- weights
  }

  # Estimate the quantile
  q <- unname(csquantile(y, weights = weights_internal, probs = probs,
                         type = type, na.rm = na.rm))

  # Estimate bandwidth using Silverman's rule for skewed distributions
  q75 <- csquantile(y, weights = weights_internal, probs = 0.75,
                    type = type, na.rm = na.rm)
  q25 <- csquantile(y, weights = weights_internal, probs = 0.25,
                    type = type, na.rm = na.rm)
  h <- 0.79 * unname(q75 - q25) * hatN^(- 1/5)

  # Estimate density at the quantile using Gaussian kernel
  fde <- (1 / (hatN * h * sqrt(2 * pi) )) *
    sum(weights_internal * exp(- (rep(q, length(y)) - y)^2 / (2 * h^2)))


  # Compute influence function for each observation
  indicator <- ifelse(y <= q, 1, 0)
  IF <- (probs - indicator) * (1 / (hatN * fde))
  IF <- unname(IF)

  return(IF)
}
