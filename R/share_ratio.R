#' Quantile-Based Share Ratio
#'
#' Estimates a quantile-based share ratio (QBSR) for measuring inequality
#' from simple or complex survey data.
#'
#' @param y A numeric vector of strictly positive values (e.g. income, wealth).
#' @param weights A numeric vector of sampling weights. If \code{NULL},
#'   all observations are equally weighted.
#' @param type Quantile estimation type: integer \code{4}--\code{9} or
#'   \code{"HD"} for Harrell--Davis (default: \code{6}). See \code{\link{csquantile}}.
#' @param prob_numerator Numeric in \eqn{(0,1)}; quantile order for the
#'   numerator (default: \code{0.80}, corresponding to the QSR top share).
#' @param prob_denominator Numeric in \eqn{(0,1)}; quantile order for the
#'   denominator (default: \code{0.20}, corresponding to the QSR bottom share).
#' @param na.rm Logical; remove missing values before computing? Default:
#'   \code{TRUE}.
#'
#' @returns A scalar numeric value representing the estimated share ratio.
#'
#' @details
#' Consider a random sample \eqn{s} of size \eqn{n}, and let \eqn{y_j} and \eqn{w_j},
#' \eqn{j \in s}, define the observed value and the sampling weight associated to the \eqn{j}-th
#' individual. Define \eqn{p_n} and \eqn{p_d} as the orders of the numerator
#' and denominator quantiles, respectively. The QBSR estimator is defined as:
#'
#' \deqn{
#'   \widehat{QBSR} =
#'   \frac{
#'     \sum_{j \in s} w_j y_j \mathbf{1}\left\{ y_j \geq \widehat{Q}(p_n) \right\}
#'   }{
#'     \sum_{j \in s} w_j y_j \mathbf{1}\left\{ y_j \leq \widehat{Q}(p_d) \right\}
#'   }
#' }
#'
#' where \eqn{\widehat{Q}(p)} is computed via \code{\link{csquantile}}, which
#' accounts for sampling weights and the specified quantile type.
#'
#' The most well-known special cases are the quintile share ratio
#' (QSR; \insertCite{langel2011quintile;textual}{inequantiles}),
#' obtained with \eqn{p_n = 0.80} and \eqn{p_d = 0.20}, and the Palma index
#' (\insertCite{palma2006globalizing;textual}{inequantiles};
#' \insertCite{palma2011homogeneous;textual}{inequantiles}),
#' obtained with \eqn{p_n = 0.90} and \eqn{p_d = 0.40}.
#'
#' @examples
#' data(synthouse)
#' eq <- synthouse$eq_income
#' w  <- synthouse$weight
#'
#' # QSR (default: top 20% vs bottom 20%)
#' share_ratio(y = eq, weights = w)
#'
#' # Palma index (top 10% vs bottom 40%)
#' share_ratio(y = eq, weights = w, prob_numerator = 0.90, prob_denominator = 0.40)
#'
#' # Compare across macro-regions (NUTS1)
#' tapply(1:nrow(synthouse), synthouse$NUTS1, function(idx) {
#'   share_ratio(y = synthouse$eq_income[idx],
#'               weights = synthouse$weight[idx])
#' })
#'
#' @references
#'
#' \insertRef{langel2011quintile}{inequantiles}
#'
#' \insertRef{palma2006globalizing}{inequantiles}
#'
#' \insertRef{palma2011homogeneous}{inequantiles}
#'
#' @seealso \code{\link{csquantile}} for quantile estimation
#'
#'
#' @family inequality indicators based on quantiles
#'
#' @export
share_ratio <- function(y, weights = NULL, type = 6, na.rm = TRUE,
                        prob_numerator   = 0.80,
                        prob_denominator = 0.20) {

  if (is.null(weights)) {
    weights <- rep(1, length(y))
  }

  if (length(y) != length(weights)) {
    stop("'y' and 'weights' must have the same length.")
  }

  # Compute quantiles
  q_n <- unname(csquantile(y, prob_numerator,   weights = weights, type = type, na.rm = na.rm))
  q_d <- unname(csquantile(y, prob_denominator, weights = weights, type = type, na.rm = na.rm))

  # Weighted income shares
  numerator   <- sum(weights * y * (y >= q_n), na.rm = na.rm)
  denominator <- sum(weights * y * (y <= q_d), na.rm = na.rm)

  if (denominator == 0) {
    warning(sprintf(
      "Denominator is zero: no observations in the bottom %.0f%%. Returning NA.",
      prob_denominator * 100
    ))
    return(NA_real_)
  }

  numerator / denominator
}
