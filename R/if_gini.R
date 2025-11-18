#' Influence Function for the Gini Coefficient
#'
#' Computes the influence function for the Gini coefficient, useful for
#' variance estimation in complex survey designs. The function implements
#' the linearization approach described in \insertCite{langel2013variance;textual}{inequantiles}.
#'
#' @param y Numeric vector of income or variable of interest.
#' @param weights Numeric vector of sampling weights. If \code{NULL} (default),
#'   equal weights are assumed (simple random sampling).
#' @param na.rm Logical. Should missing values be removed? Default is \code{TRUE}.
#'
#' @return A numeric vector of the same length as \code{y} containing the
#'   influence function values for each observation.
#'
#' @details
#' The influence function for the Gini coefficient is computed using the
#' linearization method, follwing \insertCite{deville1999variance;textual}{inequantiles}
#' framework and as defined by \insertCite{langel2013variance;textual}{inequantiles}.
#' For observation \eqn{k} with value \eqn{y_k} and
#' weight \eqn{w_k}, the influence function is:
#'
#' \deqn{z_k = \frac{2W_k(y_k - \bar{Y}_k) + \hat{Y} - \hat{N}y_k - G(\hat{Y} + y_k\hat{N})}{\hat{N}\hat{Y}}}
#'
#' where:
#' \itemize{
#'   \item \eqn{W_k = \sum_{i=1}^k w_i} is the cumulative sum of weights up to rank \eqn{k}
#'   \item \eqn{\bar{Y}_k =\frac{\sum_{l \in S} w_l y_l 1\left(W_l \leqslant W_j\right)}{W_k}} is the weighted mean of values up to rank \eqn{k}
#'   \item \eqn{\hat{N} = \sum_i w_i} is the total sum of weights
#'   \item \eqn{\hat{Y} = \sum_i w_i y_i} is the weighted total of the variable
#'   \item \eqn{G} is the Gini coefficient estimate
#' }
#'
#' The observations are ranked by their values before computing the influence function.
#'
#' @references
#'
#' \insertRef{deville1999variance}{inequantiles}
#'
#' \insertRef{langel2013variance}{inequantiles}
#'
#' @seealso \code{\link{if_qsr}} for the quintile share ratio influence function
#' and \code{\link{if_qri}} for the quantile ratio index influence function.
#'
#' @examples
#'
#' data(synthouse)
#'
#' eq <- synthouse$eq_income ### Income data
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

  # Store original order
  n <- length(y)
  original_order <- seq_len(n)

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
    G <- laeken::gini(inc = y)$value / 100

    # Influence function (equation 12 from Langel & Tillé 2013)
    z <- (2 * Nk * (x - Yk_mu) + Y - N * x - G * (Y + x * N)) / (N * Y)

    # Return in original order
    output <- z[order(r)]

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
    G <- laeken::gini(inc = y, weights = weights)$value / 100

    # Influence function (equation 12 from Langel & Tillé 2013)
    numerator <- 2 * Nk * (x - Yk_mu) + hatY - hatN * x -
      G * (hatY + x * hatN)
    z <- numerator / (hatN * hatY)

    # Return in original order
    output <- z[order(r)]
  }

  return(output)
}
