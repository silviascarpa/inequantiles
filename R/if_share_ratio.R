#' Influence Function for Quantile-Based Share Ratios
#'
#' Computes the linearized variable (influence function) for the
#' quantile-based share ratio (QBSR) using the linearization approach of
#' \insertCite{deville1999variance;textual}{inequantiles} and the derivation
#' in \insertCite{langel2011quintile;textual}{inequantiles}.
#'
#' @param y A numeric vector of strictly positive values (e.g. income, wealth).
#' @param weights A numeric vector of sampling weights. If \code{NULL},
#'   all observations are equally weighted.
#' @param type Quantile estimation type: integer 4--9 or \code{"HD"} for
#'   Harrell--Davis (default: \code{6}). See \code{\link{csquantile}}.
#' @param prob_numerator Numeric in \eqn{(0,1)}; quantile order for the
#'   numerator (default: \code{0.80}).
#' @param prob_denominator Numeric in \eqn{(0,1)}; quantile order for the
#'   denominator (default: \code{0.20}).
#' @param na.rm Logical; remove missing values before computing? Default:
#'   \code{TRUE}.
#'
#' @returns A numeric vector of the same length as \code{y} containing the
#'   linearized variable \eqn{\widehat{z}_k} for each observation.
#'
#' @details
#' \insertCite{langel2011quintile;textual}{inequantiles} derived the influence
#' function for the quintile share ratio, which generalises to any QBSR.
#' Define \eqn{p_n} and \eqn{p_d} as the quantile orders for the numerator
#' and denominator, respectively. The linearized variable is:
#'
#' \deqn{
#'   \widehat{z}_k =
#'   \frac{y_k - \widehat{I}(\widehat{Y}_{p_n})_k}{\widehat{Y}_{p_d}}
#'   -
#'   \frac{(\widehat{Y} - \widehat{Y}_{p_n})\,
#'         \widehat{I}(\widehat{Y}_{p_d})_k}{\widehat{Y}_{p_d}^2}
#' }
#'
#' where the influence function of the partial total
#' \eqn{\widehat{Y}_p = \sum_{j \in s} w_j y_j \mathbf{1}[y_j \leq \widehat{Q}(p)]}
#' is:
#'
#' \deqn{
#'   \widehat{I}(\widehat{Y}_p)_k =
#'   p\,\widehat{Q}(p) - \bigl(\widehat{Q}(p) - y_k\bigr)
#'   \mathbf{1}\bigl[y_k \leq \widehat{Q}(p)\bigr]
#' }
#'
#' and \eqn{\widehat{Y} = \sum_{j \in s} w_j y_j} is the estimated total.
#'
#' @examples
#' data(synthouse)
#' eq <- synthouse$eq_income
#' w  <- synthouse$weight
#'
#' # QSR influence function (default: p_n = 0.80, p_d = 0.20)
#' z <- if_share_ratio(eq, weights = w)
#'
#' # Palma influence function (p_n = 0.90, p_d = 0.40)
#' z_palma <- if_share_ratio(eq, weights = w,
#'                            prob_numerator = 0.90, prob_denominator = 0.40)
#'
#' @references
#' \insertRef{deville1999variance}{inequantiles}
#' \insertRef{langel2011quintile}{inequantiles}
#'
#' @seealso \code{\link{share_ratio}}, \code{\link{csquantile}}
#'
#' @family influence functions
#'
#' @export
if_share_ratio <- function(y,
                           weights          = NULL,
                           type             = 6,
                           prob_numerator   = 0.80,
                           prob_denominator = 0.20,
                           na.rm            = TRUE) {

  if (na.rm) {
    keep <- !is.na(y)
    if (!is.null(weights)) keep <- keep & !is.na(weights)
    y <- y[keep]
    if (!is.null(weights)) weights <- weights[keep]
  }

  if (is.null(weights)) {
    weights <- rep(1, length(y))
  }

  if (length(y) != length(weights)) {
    stop("'y' and 'weights' must have the same length.")
  }

  # Weighted quantiles
  Q_n <- unname(csquantile(y, prob_numerator,   weights = weights, type = type, na.rm = FALSE))
  Q_d <- unname(csquantile(y, prob_denominator, weights = weights, type = type, na.rm = FALSE))

  # Weighted partial totals and grand total
  Y_n   <- sum(weights * y * (y <= Q_n))
  Y_d   <- sum(weights * y * (y <= Q_d))
  Y_tot <- sum(weights * y)

  # Influence function of partial total Y_p (vectorised over all observations)
  # I(Y_p)_k = p * Q(p) - (Q(p) - y_k) * 1[y_k <= Q(p)]
  I_Yn <- prob_numerator   * Q_n - (Q_n - y) * (y <= Q_n)
  I_Yd <- prob_denominator * Q_d - (Q_d - y) * (y <= Q_d)

  # Delta method: IF of (Y_top / Y_d) where Y_top = Y_tot - Y_n
  z_k <- (y - I_Yn) / Y_d - (Y_tot - Y_n) * I_Yd / Y_d^2

  unname(z_k)
}
