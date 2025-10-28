#' Estimator of Quantile in case of Simple and Complex Sampling Design
#'
#' Computes quantiles for weighted or unweighted data, allowing for sampling weights
#' and several interpolation types.
#' The method extends the standard quantile
#' definitions of \insertCite{hyndman1996sample;textual}{inequantiles} and
#' \insertCite{harrell1982new;textual}{inequantiles} estimator to the case of complex survey data
#' by incorporating sampling weights into the cumulative distribution function (CDF)
#' and interpolation points, as proposed in \insertCite{scarpa2025inference;textual}{inequantiles}.
#'
#' @param y Numeric vector of observations.
#' @param weights Optional numeric vector of sampling weights (default: NULL for equal weights)
#' @param probs Numeric vector of probabilities (default = seq(0, 1, 0.1)).
#' @param type Quantile estimation algorithm: integer 4-9 or "HD" for Harrell-Davis (default: 4)
#' @param na.rm Logical indicating whether to remove NA values (default: TRUE)
#'
#' @details
#' Consider a random sample \eqn{s} of size \eqn{n}. Let \eqn{Y_1, \ldots, Y_n} be the sample observations from a finite population,
#' with order statistics \eqn{Y_{(1)} \le \ldots \le Y_{(n)}} and corresponding sampling
#' weights \eqn{w_1, \ldots, w_n}. Define the cumulative weights
#' \eqn{W_j = \sum_{i \le j} w_i} and the total weight \eqn{W_n = \sum_{i=1}^n w_i}.
#' The weighted quantile estimator is computed as a linear interpolation between
#' adjacent order statistics:
#'
#' \deqn{
#'   \widehat{Q}(p)
#'   = Y_{(k-1)} +
#'     (Y_{(k)} - Y_{(k-1)})
#'     \frac{p - \widehat{r}_{k-1}}{\widehat{r}_k - \widehat{r}_{k-1}},
#' }
#'
#' where \eqn{\widehat{r}_k} denotes the estimated cumulative distribution function
#' (the “plotting position”), and the order \eqn{k} is such that
#' \eqn{W_{k-1} - m_{k-1} < W_n p < W_k - m_k},
#' with \eqn{m_k} determined by the interpolation method.
#'
#' The function supports several interpolation rules (types 4–9), extending the
#' quantile definitions in \insertCite{hyndman1996sample;textual}{inequantiles} to incorporate sampling weights.
#'
#' The table below summarizes the six interpolation types (4–9) extended from
#' \insertCite{hyndman1996sample;textual}{inequantiles} to incorporate sampling weights,
#' as described in \insertCite{scarpa2025inference;textual}{inequantiles}.
#'
#' \tabular{llll}{
#' \strong{Type} \tab \strong{Estimator \eqn{\widehat{r}_k}} \tab
#' \strong{Interpolation \eqn{\widehat{m}_k}} \tab
#' \strong{Selection rule for \eqn{k}} \cr
#' 4 \tab \eqn{W_k / W_n} \tab 0 \tab \eqn{W_{k-1} \le W_n p < W_k} \cr
#' 5 \tab \eqn{(W_k - \tfrac{1}{2} w_k) / W_n} \tab \eqn{w_k / 2} \tab
#'   \eqn{W_{k-1} - \tfrac{w_{k-1}}{2} \le W_n p < W_k - \tfrac{w_k}{2}} \cr
#' 6 \tab \eqn{W_k / (W_n + w_n)} \tab \eqn{w_n p} \tab
#'   \eqn{W_{k-1} \le (W_n + w_n)p < W_k} \cr
#' 7 \tab \eqn{W_{k-1} / W_{n-1}} \tab \eqn{w_k - w_n p} \tab
#'   \eqn{W_{k-2} \le W_{n-1}p < W_{k-1}} \cr
#' 8 \tab \eqn{(W_k - \tfrac{1}{3}w_k) / (W_n + \tfrac{w_n}{3})} \tab
#'   \eqn{\tfrac{w_k}{3} + \tfrac{w_n}{3}p} \tab
#'   \eqn{W_{k-1} - \tfrac{w_{k-1}}{3} \le (W_n - \tfrac{w_n}{3})p < W_k - \tfrac{w_k}{3}} \cr
#' 9 \tab \eqn{(W_k - \tfrac{3}{8}w_k) / (W_n + \tfrac{1}{4}w_n)} \tab
#'   \eqn{\tfrac{3}{8}w_k + \tfrac{w_n}{4}p} \tab
#'   \eqn{W_{k-1} - \tfrac{3w_{k-1}}{8} \le (W_n + \tfrac{w_n}{4})p < W_k - \tfrac{3w_k}{8}} \cr
#' }
#'
#' For unweighted data, the function returns the standard R quantiles.
#'
#' The Harrell–Davis estimator ("HD") is also extended to the weighted case as
#' proposed in Kreutzmann (2018), by redefining the weighting coefficients
#' \eqn{\widehat{\mathcal{W}}_j(p)} for order statistics as:
#'
#' \deqn{
#' \widehat{\mathcal{W}}_j(p)
#'   = b_{(W_j / W_n)}\{(W_n + w_n)p, W_n - (W_n + w_n)p + w_n\}
#'   - b_{(W_{j-1}/W_n)}\{(W_n + w_n)p, W_n - (W_n + w_n)p + w_n\},
#' }
#'
#' where \eqn{b_x(a,b)} denotes the incomplete beta function.
#'
#' The resulting quantile estimator is
#' \eqn{\widehat{Q}_{HD}(p) = \sum_{j \in s} \widehat{\mathcal{W}}_j(p) Y_{(j)}.}
#' For unweighted data, the function returns the Harrell & Davis quantile estimator.
#'
#'
#' @examples
#'
#' data(synthouse)
#'
#' y <- synthouse$eq_income ### Income data
#'
#' # Compute unweighted quantiles (default: type = 4)
#' csquantile(y, probs = c(0.25, 0.5, 0.75), type = 6)
#'
#' # Consider the sampling weights
#' w <- synthouse$weight
#' csquantile(y, weight = w, probs = c(0.25, 0.5, 0.75), type = 6)
#'
#'
#'
#' @return
#' A named numeric vector of estimated quantiles corresponding to \code{probs}.
#'
#' @references
#' \insertRef{harrell1982new}{inequantiles}
#'
#' \insertRef{hyndman1996sample}{inequantiles}
#'
#' \insertRef{scarpa2025inference}{inequantiles}
#'
#' @importFrom stats integrate pbeta quantile var
#'
#' @export

csquantile <- function(y,
                       weights = NULL,
                       probs = seq(0, 1, 0.1),
                       type = 4,
                       na.rm = FALSE) {

  # Handle missing values
  if (na.rm) {
    na_idy <- is.na(y)
    if (any(na_idy)) {
      y <- y[!na_idy]
      if (!is.null(weights)) {
        weights <- weights[!na_idy]
      }
    }
  }

  # Sort data and weights
  order_idy <- order(y)
  y <- y[order_idy]
  n <- length(y)

  # ============================================================================
  # HARRELL-DAVIS ESTIMATOR
  # ============================================================================
  if (identical(type, "HD")) {

    if (is.null(weights)) {
      # Unweighted HD
      m <- n + 1
      ps <- probs[probs > 0 & probs < 1]

      a <- outer((0:n) / n, ps,
                 function(y, p, m) pbeta(y, p * m, (1 - p) * m),
                 m = m)

    } else {
      # Weighted HD
      w <- weights[order_idy]
      total_weight <- sum(w)

      if (total_weight < 2) {
        return(rep(NA, length(probs)))
      }

      m <- total_weight + 1
      ps <- probs[probs > 0 & probs < 1]
      cum_prop <- c(0, cumsum(w)) / total_weight

      a <- outer(cum_prop, ps,
                 function(y, p, m) pbeta(y, p * m, (1 - p) * m),
                 m = m)
    }

    # Compute weights and quantiles
    w_matriy <- a[-1, , drop = FALSE] - a[-(n + 1), , drop = FALSE]
    r <- drop(y %*% w_matriy)

    # Handle boundary probabilities
    rp <- range(probs)
    pp <- ps

    if (rp[1] == 0) {
      r <- c(y[1], r)
      pp <- c(0, pp)
    }

    if (rp[2] == 1) {
      r <- c(r, y[n])
      pp <- c(pp, 1)
    }

    r <- r[match(probs, pp)]
    names(r) <- paste0(probs * 100, "%")

    return(r)
  }

  # ============================================================================
  # STANDARD QUANTILE TYPES (4-9)
  # ============================================================================

  # Unweighted case: use standard R quantile
  if (is.null(weights)) {
    q <- quantile(y, type = type, probs = probs, na.rm = FALSE)
    return(q)
  }

  # Weighted case
  w <- weights[order_idy]
  total_weight <- sum(w)
  cumw <- cumsum(w)

  q <- sapply(probs, function(p) {

    # Boundary cases
    if (p == 0) return(y[1])
    if (p > 1 - 4 * .Machine$double.eps) return(y[n])

    # Type-specific computations for l and interpolation parameters
    if (type == 4) {
      l <- max(which(cumw <= total_weight * p), 2)
      u <- min(n, l + 1)
      condition <- total_weight * p > cumw[l]
      numerator <- total_weight * p - cumw[l]
      denominator <- w[u]

    } else if (type == 5) {
      l <- max(2, which(cumw <= w / 2 + total_weight * p))
      if (l == length(w)) return(y[n])
      u <- min(n, l + 1)
      condition <- total_weight * p > (cumw[l] - w[l] / 2)
      numerator <- total_weight * p - (cumw[l] - w[l] / 2)
      denominator <- (w[u] + w[l]) / 2

    } else if (type == 6) {
      adjusted_weight <- total_weight + w[n]
      l <- max(2, which(cumw <= adjusted_weight * p))
      u <- min(n, l + 1)
      condition <- adjusted_weight * p > cumw[l]
      numerator <- adjusted_weight * p - cumw[l]
      denominator <- w[u]

    } else if (type == 7) {
      adjusted_weight <- total_weight - w[n]
      l <- max(2, max(which(cumw <= w + adjusted_weight * p)))
      u <- min(n, l + 1)
      condition <- adjusted_weight * p > cumw[l - 1]
      numerator <- adjusted_weight * p - cumw[l - 1]
      denominator <- w[l]

    } else if (type == 8) {
      adjusted_weight <- total_weight + w[n] / 3
      l <- max(2, which(cumw <= w / 3 + adjusted_weight * p))
      u <- min(n, l + 1)
      condition <- adjusted_weight * p > (cumw[l] - w[l] / 3)
      numerator <- adjusted_weight * p - (cumw[l] - w[l] / 3)
      denominator <- 2 * w[u] / 3 + w[l] / 3

    } else if (type == 9) {
      adjusted_weight <- total_weight + w[n] / 4
      l <- max(2, which(cumw <= 3 * w / 8 + adjusted_weight * p))
      u <- min(n, l + 1)
      condition <- adjusted_weight * p > (cumw[l] - 3 * w[l] / 8)
      numerator <- adjusted_weight * p - (cumw[l] - 3 * w[l] / 8)
      denominator <- 5 * w[u] / 8 + 3 * w[l] / 8

    } else {
      stop("Unsupported quantile type. Use types 4-9 or 'HD'.")
    }

    # Common computation for all types
    p_ok <- !is.na(p)
    i <- which(!p_ok | (condition & y[l] != y[u]))
    h <- if (length(i) > 0) numerator / denominator else 0

    return((1 - h) * y[l] + h * y[u])
  })

  names(q) <- paste0(probs * 100, "%")
  return(q)
}
