#' Influence Function for the quantile ratio index
#'
#' Computes the influence function of the QRI for all observations.
#'
#' @param y A numeric vector of data values
#' @param weights A numeric vector of sampling weights (optional)
#' @param k Optional integer, index of the unit for which to compute the influence function.
#'   If NULL (default), the function returns the IF for all observations.
#' @param type Quantile estimation type: integer 4–9 or "HD" for Harrell–Davis (default: 7)
#'
#' @return A numeric value (if k is provided) or a numeric vector (if k is NULL) of influence function values
#'
#' @details
#' The influence function is computed as:
#' \deqn{I(QRI)_i = - \int_0^1 \frac{IF(Q(p/2))_i Q(1-p/2) - IF(Q(1-p/2))_i Q(p/2)}{Q(1-p/2)^2} dp}
#'
#' The density is estimated using a Gaussian kernel with bandwidth
#' \eqn{h = 0.79 \cdot IQR \cdot \hat{N}^{-1/5}}.
#'
#' #' #' @references
#' Scarpa, S., Ferrante, M.R., & Sperlich, S. (2025). Inference for the Quantile Ratio
#'   Inequality Index in the Context of Survey Data. \emph{Journal of Survey Statistics and Methodology}
#'
#' @export

if_qri <- function(y, weights = NULL, type = 6) {

  n <- length(y)

  # Determine effective sample size
  if (is.null(weights)) {
    hatN <- n
    w <- rep(1, n)
  } else {
    hatN <- sum(weights)
    w <- weights
  }

  # Compute bandwidth using IQR
  q75 <- csquantile(y, weights = weights, probs = 0.75, type = type, na.rm = FALSE)
  q25 <- csquantile(y, weights = weights, probs = 0.25, type = type, na.rm = FALSE)
  h <- 0.79 * unname(q75 - q25) * hatN^(-1/5)

  # Density estimator function
  fde <- function(r) {
    density <- (1 / (hatN * h * sqrt(2 * pi))) *
      sum(w * exp(-(rep(r, n) - y)^2 / (2 * h^2)))
    return(density)
  }

  # Ratio influence function for observation k
  ratio_if_k <- function(prob, k) {
    # Compute quantiles
    qlow <- unname(csquantile(y, weights = weights, probs = prob / 2,
                              type = type, na.rm = FALSE))
    qhigh <- unname(csquantile(y, weights = weights, probs = 1 - prob / 2,
                               type = type, na.rm = FALSE))

    # Compute density at quantiles
    f_low <- fde(qlow)
    f_high <- fde(qhigh)

    # Influence function for lower quantile
    IF_low <- - (1 / (hatN * f_low)) * (as.numeric(y[k] <= qlow) - prob / 2)

    # Influence function for upper quantile
    IF_high <- - (1 / (hatN * f_high)) * (as.numeric(y[k] <= qhigh) - (1 - prob / 2))

    # Ratio influence function
    ratio_if_value <- (IF_low * qhigh - IF_high * qlow) / (qhigh^2)

    return(ratio_if_value)
  }

  # Compute IF for all observations
  IF_values <- sapply(1:n, function(k) {
    result <- integrate(function(p) ratio_if_k(p, k), lower = 0, upper = 1)
    return(- result$value)
  })

  return(IF_values)
}
