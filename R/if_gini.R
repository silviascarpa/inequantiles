# Internal helper: weighted Gini coefficient.
# Formula: \insertCite{langel2013variance;textual}{inequantiles}, equation 6.
# Not exported; used by if_gini() and inequantiles().
.gini_coef <- function(y, weights = NULL) {
  if (is.null(weights)) weights <- rep(1L, length(y))
  ord  <- order(y)
  y    <- y[ord]
  w    <- weights[ord]
  hatN <- sum(w)
  hatY <- sum(w * y)
  Wj   <- cumsum(w)
  (2 * sum(w * y * Wj) - sum(w^2 * y)) / (hatN * hatY) - 1
}


#' Influence Function for the Gini Coefficient
#'
#' Computes the influence function for the Gini coefficient, useful for
#' variance estimation and linearization in complex survey designs \insertCite{langel2013variance;textual}{inequantiles}.
#'
#' @param y Numeric vector of income or variable of interest.
#' @param weights Numeric vector of sampling weights. If \code{NULL} (default),
#'   equal weights are assumed (simple random sampling).
#' @param na.rm Logical. Should missing values be removed? Default is \code{TRUE}.
#'
#' @returns A numeric vector of the same length as \code{y} containing the
#'   influence function values for each observation, returned in the same order
#'   as the input \code{y}.
#'
#' @details
#' The influence function for the Gini coefficient is computed using the
#' linearization method, following \insertCite{deville1999variance;textual}{inequantiles}
#' framework and as defined by \insertCite{langel2013variance;textual}{inequantiles}.
#' The influence function for Gini is:
#'
#' \deqn{{I}(\widehat{G})_{k} = \frac{2W_k(y_k - \bar{Y}_k) + \hat{Y} - \hat{N}y_k - G(\hat{Y} + y_k\hat{N})}{\hat{N}\hat{Y}}}
#'
#' where:
#' \itemize{
#'   \item \eqn{W_k = \sum_{i=1}^k w_i} is the cumulative sum of weights up to rank \eqn{k}
#'   \item \eqn{\bar{Y}_k =\frac{\sum_{l \in S} w_l y_l 1\left(W_l \leqslant W_k\right)}{W_k}} is the weighted mean of values up to rank \eqn{k}
#'   \item \eqn{\hat{N} = \sum_i w_i} is the total sum of weights
#'   \item \eqn{\hat{Y} = \sum_i w_i y_i} is the weighted total of the variable
#'   \item \eqn{G} is the Gini coefficient estimate
#' }
#'
#' @references
#'
#' \insertRef{deville1999variance}{inequantiles}
#'
#' \insertRef{langel2013variance}{inequantiles}
#'
#'
#' @family influence functions
#'
#' @examples
#'
#' data(synthouse)
#'
#' eq <- synthouse$eq_income # Equivalized disposable income
#'
#' # Simple example
#' z <- if_gini(eq)
#'
#' # With weights
#' w <- synthouse$weight
#' z_weighted <- if_gini(y = eq, weights = w)
#'
#'
#'
#' @export
if_gini <- function(y, weights = NULL, na.rm = TRUE) {

  # Handle missing values
  if (na.rm) {
    if (is.null(weights)) {
      valid <- !is.na(y)
      y <- y[valid]
    } else {
      valid <- !is.na(y) & !is.na(weights)
      y <- y[valid]
      weights <- weights[valid]
    }
  }

  # Check for any remaining NAs
  if (any(is.na(y))) {
    stop("Missing values in y. Set na.rm = TRUE to remove them.")
  }
  if (!is.null(weights) && any(is.na(weights))) {
    stop("Missing values in weights. Set na.rm = TRUE to remove them.")
  }

  n <- length(y)

  # Case 1: Unweighted (simple random sampling)
  if (is.null(weights)) {

    # Rank and sort
    r <- rank(y, ties.method = "first")
    ord <- order(y)
    x <- y[ord]

    N <- length(x)
    Y <- sum(x)
    Nk <- seq_len(N)
    Yk <- cumsum(x)
    Yk_mu <- Yk / Nk

    # Gini coefficient
    G <- .gini_coef(y)

    # Influence function (equation 12 from Langel & Tillé 2013)
    # z[k] = IF for the k-th smallest observation
    z <- (2 * Nk * (x - Yk_mu) + Y - N * x - G * (Y + x * N)) / (N * Y)

    # Return in original order: obs i has rank r[i], so its IF is z[r[i]]
    output <- z[r]

  } else {
    # Case 2: Weighted

    # Check weights
    if (length(weights) != length(y)) {
      stop("Length of weights must equal length of y.")
    }
    if (any(weights < 0)) {
      stop("Weights must be non-negative.")
    }

    # Rank and sort
    r <- rank(y, ties.method = "first")
    ord <- order(y)
    x <- y[ord]
    w <- weights[ord]

    # Weighted totals
    hatN <- sum(w)
    hatY <- sum(w * x)

    # Cumulative sums
    Nk <- cumsum(w)
    Yk <- cumsum(w * x)
    Yk_mu <- Yk / Nk

    # Gini coefficient
    G <- .gini_coef(y, weights)

    # Influence function (equation 12 from Langel & Tillé 2013)
    # z[k] = IF for the k-th smallest observation
    numerator <- 2 * Nk * (x - Yk_mu) + hatY - hatN * x -
      G * (hatY + x * hatN)
    z <- numerator / (hatN * hatY)

    # Return in original order: obs i has rank r[i], so its IF is z[r[i]]
    output <- z[r]
  }

  return(output)
}
