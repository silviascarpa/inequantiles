#' Influence Function for Quantiles
#'
#' Computes the influence function for quantile estimators under simple and complex survey designs
#'
#' @param y A numeric vector of data values
#' @param weights A numeric vector of sampling weights (optional)
#' @param probs A numeric value specifying the probability for the quantile (e.g., 0.5 for median)
#' @param type Quantile estimation type: integer 4-9 or "HD" for Harrell-Davis (default: 6)
#' @param na.rm Logical, should missing values be removed? (default: TRUE)
#'
#' @return A numeric vector containing the influence function values for each observation
#'
#' @details
#' The influence function of a quantile Q(p) is given by:
#' \deqn{I(Q(p))_i = \frac{p - \mathbbm{1}(y_i \leq Q(p))}{F'(Q(p)) N}}
#'
#' where F'(Q(p)) is estimated using a Gaussian kernel density estimator with bandwidth
#' h = 0.79 * IQR * N^{-1/5} as suggested by Verma (2011).
#'
#' #' @references
#' Osier, G., (2009), “Variance estimation for complex indicators of poverty and inequality using
#'  linearization techniques”, Survey Research Methods, 3, 167–195
#' Scarpa, S., Ferrante, M.R., & Sperlich, S. (2025). Inference for the Quantile Ratio
#'   Inequality Index in the Context of Survey Data. \emph{Journal of Survey Statistics and Methodology}.
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

  return(IF)
}
