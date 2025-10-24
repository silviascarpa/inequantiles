#' Influence Function for the quantile ratio index
#'
#' Computes the influence function of the quantile ratio index (QRI) in the context
#' of finite population (see Deville, 1999, for a definition) for all observations,
#' as defined in Scarpa et al. (2025), under simple and complex sampling.
#'
#' @param y A numeric vector of data values
#' @param weights A numeric vector of sampling weights (optional). If \code{NULL}, equal weights are assumed.
#' @param type Quantile estimation type: integer 4–9 or HD for Harrell–Davis (default: 6)
#'         See \code{csquantile} documentation for a complete description.
#' @return A numeric vector of influence function values (one for each observation)
#'
#' @details
#' The influence function for the QRI is computed on each observation as
#' \deqn{
#'   \widehat{z}_i = - \int_0^1
#'   \frac{
#'   \left(
#'     \frac{\frac{p}{2} - \mathbf{1}(y_i \leq \widehat{Q}(p/2))}
#'     {\widehat{f}(\widehat{Q}(p/2)) \, \widehat{N}}
#'   \right)
#'   \widehat{Q}(1 - p/2)
#'   -
#'   \left(
#'     \frac{(1 - \frac{p}{2}) - \mathbf{1}(y_i \leq \widehat{Q}(1 - p/2))}
#'     {\widehat{f}(\widehat{Q}(1 - p/2)) \, \widehat{N}}
#'   \right)
#'   \widehat{Q}(p/2)
#'   }{
#'     \widehat{Q}(1 - p/2)^2
#'   } \, dp
#' }
#'
#' where:
#' \itemize{
#'   \item \eqn{\widehat{Q}(p)} is the weighted sample quantile of order \eqn{p},
#'     computed using the internal function \code{csquantile()},
#'   \item \eqn{\widehat{f}(\cdot)} denotes the estimated income density function,
#'   \item \eqn{\widehat{N} = \sum_i w_i} is the estimated population size.
#' }
#'
#' The density function \eqn{\widehat{f}(y)} is estimated via a Gaussian kernel smoother:
#' \deqn{
#'   \widehat{f}(y) =
#'   \frac{1}{\widehat{N}} \sum_{j \in s} w_j
#'   K\!\left(\frac{y - y_j}{h}\right)
#'   =
#'   \frac{1}{\widehat{N}\, h \sqrt{2\pi}}
#'   \sum_{j \in s} w_j
#'   \exp\!\left\{ -\frac{(y - y_j)^2}{2h^2} \right\},
#' }
#' where \eqn{K(\cdot)} is the Gaussian kernel.
#'
#' The bandwidth is chosen as:
#' \deqn{
#'   h = 0.79 \cdot \mathrm{IQR} \cdot \widehat{N}^{-1/5},
#' }
#' where \eqn{\mathrm{IQR}} is the interquartile range of the weighted sample.
#'
#' @examples
#'
#' # On synthetic data
#' eq_synth <- rlnorm(30, 9, 0.7)
#' IF_synth <- if_qri(y = eq_synth)
#'
#' # On real data
#' eq <- synthouse$eq_income[1:30] ## Take some observations (as example)
#' w <- synthouse$weight[1:30]
#' IF_qri <- if_qri(y = eq, weights = w, type = 6)
#'
#'
#'
#' @references
#'
#' Deville, J.C., (1999), “Variance estimation for complex statistics and estimators:
#' Linearization and residual techniques”, *Survey methodology*, 25, 193–204
#'
#' Scarpa, S., Ferrante, M.R., & Sperlich, S. (2025). "Inference for the Quantile Ratio
#'   Inequality Index in the Context of Survey Data".
#'   *{Journal of Survey Statistics and Methodology}*, smaf024
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
