#' Estimator of quantile in case of simple and complex sampling design
#'
#' Computes quantiles for weighted or unweighted data, allowing several interpolation, as studied by
#' Scarpa et al (2025). Supports standard quantile types (4-9) and the
#' Harrell-Davis estimator ("HD") when no weights are considered.
#'
#' @param y Numeric vector of observations.
#' @param weights Optional numeric vector of sampling weights (default: NULL for equal weights)
#' @param probs Numeric vector of probabilities (default = seq(0, 1, 0.1)).
#' @param type Quantile algorithm: integer 4-9 or "HD" for Harrell-Davis (default: 5)
#' @param na.rm Logical indicating whether to remove NA values (default: FALSE)
#'
#'
#' @return
#' A named numeric vector of quantiles corresponding to \code{probs}.
#'
#' @references
#' Hyndman, R.J. & Fan, Y. (1996). Sample quantiles in statistical packages.
#'   \emph{The American Statistician}, 50(4), 361–365.
#'
#' Scarpa, S., Ferrante, M.R., & Sperlich, S. (2025). Inference for the Quantile Ratio
#'   Inequality Index in the Context of Survey Data. \emph{Journal of Survey Statistics and Methodology}.
#'
#'
#' @export

csquantile <- function(y,
                       weights = NULL,
                       probs = seq(0, 1, 0.1),
                       type = 5,
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
