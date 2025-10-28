#' Quantile Ratio Index Estimator in Superpopulation
#'
#' Computes the theoretical quantile ratio index (QRI) for measuring inequality
#' for a given parametric distribution.
#'
#' @param qfunction A quantile function (e.g., \code{qnorm}, \code{qlnorm}, \code{qgamma}).
#' @param lower Lower bound of integration. Default is 0.
#' @param upper Upper bound of integration. Default is 1.
#' @param subdivisions Maximum number of subintervals for integration. Default is 1000L.
#' @param ... Additional parameters to pass to \code{qfunction} (e.g., distribution parameters).
#'
#'
#'@return A numeric value representing the theoretical QRI for the specified
#'   parametric distribution. Values range from 0 (perfect equality) to 1
#'   (maximum inequality).#'
#'
#' @details
#' The QRI was proposed by \insertCite{prendergast2018simple}{inequantiles} fot the
#' economic inequality measurement. It is calculated as:
#' \deqn{QRI = 1 - \int_0^1 R(p) dp}
#' where \eqn{R(p) = Q(p/2) / Q(1 - p/2)} is the ratio between symmetric quantiles
#' function.
#'
#' This function computes the (superpopulation) QRI for
#' theoretical parametric distributions, as opposed to \code{\link{qri}} which estimates
#' the QRI from sample data.
#'
#' @examples
#' # Log-normal distribution
#' superpop_qri(qlnorm, meanlog = 9, sdlog = 0.3)
#' superpop_qri(qlnorm, meanlog = 9, sdlog = 1.4)
#'
#' # Weibull distribution
#' superpop_qri(qweibull, shape = 1.7, scale = 30000)
#' superpop_qri(qweibull, shape = 1.7, scale = 30000)
#'
#'
#' @seealso \code{\link{qri}} for the sample-based QRI estimator
#'
#' @references
#' \insertRef{prendergast2018simple}{inequantiles}
#'
#'
#' @importFrom Rdpack reprompt
#'
#'
#' @export
superpop_qri <- function(qfunction, lower = 0, upper = 1, subdivisions = 1000L, ...) {
  R <- function(p, ...) {
    qfunction(p/2, ...) / qfunction(1 - p/2, ...)
  }

  1 - integrate(R, lower = lower, upper = upper,
                subdivisions = subdivisions, ...)$value
}
