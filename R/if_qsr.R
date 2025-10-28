#' Influence Function for the Quintile Share Ratio (QSR)
#'
#' @description
#' Computes the linearized variable (influence function) for the quintile share
#' ratio (QSR) using the linearization approach based on
#' \insertCite{deville1999variance;textual}{inequantiles}
#'  and the derivation in \insertCite{langel2011quintile;textual}{inequantiles}
#'
#' @param y A numeric vector of income or other continuous variable.
#' @param weights A numeric vector of sampling weights. If \code{NULL} (default),
#'   equal weights are assumed for all observations.
#' @param type Quantile estimation type: integer 4–9 or "HD" for Harrell–Davis (default: 4)
#'
#'
#' @return A numeric vector of the same length as \code{y} containing the
#'   linearized variable \eqn{\widehat{z}_k} for each observation.
#'
#'
#' @details
#' According to the definition by \insertCite{langel2011quintile;textual}{inequantiles} ,
#' the influence function of QSR is given by:
#'
#' \deqn{I(\mathrm{QSR})_k=\frac{y_k-I\left(Y_{0.8}\right)_k}{Y_{0.2}}-\frac{\left(Y-Y_{0.8}\right) I\left(Y_{0.2}\right)_k}{Y_{0.2}^2}}
#'
#' where \eqn{I(Y_p)_k = p Q(p) - (Q(p) - y_k) \mathbb{1}[y_k \leq Q(p)]} is the
#' influence function of the partial total up to quantile \eqn{p} and
#' \eqn{Y_{\alpha} = \sum_{j} y_j\mathbb{1}[y_k \leq {Q}(\alpha)]}.
#'
#'
#' The estimated linearized variable is:
#' \deqn{\widehat{z}_k=\frac{y_k-\left\{0.8 \widehat{Q}(0.8)-\left(\widehat{Q}(0.8)-y_k\right) \mathbb{1}\left[y_k \leq \widehat{Q}(0.8)\right]\right\}}{\widehat{Y}_{0.2}}-\frac{\left(\widehat{Y}-\widehat{Y}_{0.8}\right)\left\{0.2 \widehat{Q}(0.2)-\left(\widehat{Q}(0.2)-y_k\right) \mathbb{1}\left[y_k \leq \widehat{Q}(0.2)\right]\right\}}{\widehat{Y}_{0.2}^2}}
#'
#' where \eqn{\hat{Y}_{\alpha} = \sum_{j \in s} w_j y_j\mathbb{1}[y_k \leq \widehat{Q}(\alpha)]} and
#' \eqn{\widehat{Q}(p)} is the weighted p-th quantile estimator, computed using the internal function \code{csquantile()}.
#'
#'
#' @references
#'
#' \insertRef{deville1999variance}{inequantiles}
#'
#' \insertRef{langel2011quintile}{inequantiles}
#'
#'
#' @seealso \code{\link{qsr}} for the QSR estimator, \code{\link{csquantile}}
#'   for weighted quantile estimation.
#'
#' @examples
#'
#' data(synthouse)
#'
#' eq <- synthouse$eq_income ### Income data
#'
#' # Simple example
#' z <- if_qsr(eq)
#'
#' # With weights
#' w <- synthouse$weight
#' z_weighted <- if_qsr(y = eq, weights = w)
#'
#' @export
#'
if_qsr <- function(y, weights = NULL, type = 4) {

  # Gestione pesi
  if (is.null(weights)) {
    weights <- rep(1, length(y))
  }

  # Validate that y and weights have the same length
  if (length(y) != length(weights)) {
    stop("'y' and 'weights' must have the same length")
  }

  n <- length(y)

  # Calcola i quantili pesati
  Q_80 <- csquantile(y, 0.8, weights = weights, type = type)
  Q_20 <- csquantile(y, 0.2, weights = weights, type = type)

  # Calcola le somme parziali pesate
  Y_02 <- sum(weights * y * (y <= Q_20))
  Y_08 <- sum(weights * y * (y <= Q_80))

  # Y = somma pesata totale
  Y_tot <- sum(weights * y)

  # Funzione di influenza del totale parziale
  # I(Y_p)_k = p * Q_p - (Q_p - y_k) * 1[y_k <= Q_p]
  if_partial_total <- function(y_k, Q_p, p) {
    p * Q_p - (Q_p - y_k) * (y_k <= Q_p)
  }

  # Calcola I(Y_0.8)_k e I(Y_0.2)_k per ogni k
  I_Y08 <- sapply(y, function(y_k) if_partial_total(y_k, Q_80, 0.8))
  I_Y02 <- sapply(y, function(y_k) if_partial_total(y_k, Q_20, 0.2))

  # Calcola la funzione di influenza del QSR
  # z_k = (y_k - I(Y_0.8)_k) / Y_0.2 - (Y - Y_0.8) * I(Y_0.2)_k / Y_0.2^2
  z_k <- (y - I_Y08) / Y_02 - ((Y_tot - Y_08) * I_Y02) / (Y_02^2)

  return(z_k)
}

