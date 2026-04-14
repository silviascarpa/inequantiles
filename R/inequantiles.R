#' Quantile-based inequality indicators
#'
#' Computes one or more quantile-based inequality indicators simultaneously among
#' QRI, QSR, Palma ratio, and percentile ratio, optionally with standard
#' errors (estimated from the same set of bootstrap replicates), ensuring full
#' comparability across indicators.
#'
#' @param y A numeric vector of strictly positive values (e.g. income, wealth, expenditure).
#' @param weights A numeric vector of sampling weights. If \code{NULL},
#'   all observations are equally weighted.
#' @param indicators Character vector specifying which indicators to compute.
#'   Use \code{"all"} (default) for all four, or any subset of
#'   \code{"qri"}, \code{"qsr"}, \code{"palma"}, \code{"ratio_quantiles"}.
#' @param se Logical; if \code{TRUE}, standard errors are estimated via the
#'   rescaled bootstrap; see \code{\link{rescaled_bootstrap}}. Requires \code{data} and \code{strata} (see below).
#' @param type Quantile estimation type (integer 4--9 or \code{"HD"} for
#'   Harrell-Davis). Default: \code{6}. See \code{\link{csquantile}}.
#' @param na.rm Logical; remove missing values before computing? Default:
#'   \code{TRUE}.
#' @param M Integer; number of quantile-ratio grid points for the QRI
#'   (default: \code{100}). Only used when the QRI is estimated; see \code{\link{qri}}.
#' @param prob_num Numeric in \eqn{(0,1)}; numerator quantile for the
#'   percentile ratio (default: \code{0.90}). Only used when percentiles ratio is estimated;
#'   see \code{\link{ratio_quantiles}}.
#' @param prob_den Numeric in \eqn{(0,1)}; denominator quantile for the
#'   percentile ratio (default: \code{0.10}). Only used when percentiles ratio is estimated;
#'   see \code{\link{ratio_quantiles}}.
#' @param se Logical; if \code{TRUE}, bootstrap standard errors are computed.
#' @param B Integer; number of bootstrap replicates (default: \code{200}).
#'   Only used when \code{se = TRUE}.
#' @param seed Integer; random seed for reproducibility. Only used when
#'   \code{se = TRUE}.
#' @param data A data frame containing the survey design variables (strata,
#'   PSU). Required when \code{se = TRUE}.
#' @param strata Character string; name of the stratification column in
#'   \code{data}. Required when \code{se = TRUE}.
#' @param psu Character string; name of the PSU column in \code{data}.
#'   Required for two-stage complex designs.
#' @param N_h Optional named numeric vector of stratum population sizes for
#'   the finite population correction. See \code{\link{rescaled_bootstrap}}.
#' @param m_h Optional vector of bootstrap sample sizes per stratum.
#'   Defaults to the Rao-Wu formula. See \code{\link{rescaled_bootstrap}}.
#'
#' @returns A list with components:
#'   \item{estimates}{Numeric vector of point estimates of inequality indicators.}
#'   \item{se}{Numeric vector of standard errors, or \code{NULL} when
#'     \code{se = FALSE}.}
#'   \item{B}{Number of bootstrap replicates used, or \code{NULL}.}
#'   \item{call}{The matched function call.}
#'
#' @details
#' All indicators are computed from the same specified \code{\link{csquantile}} type.
#' When \code{se = TRUE}, a \emph{single} bootstrap loop is run through rescaled
#' bootstrap method (see details in \code{\link{rescaled_bootstrap}}) and all
#' indicators are evaluated on each replicate, so standard errors are based on
#' identical resamples and are directly comparable.
#'
#'
#' @seealso \code{\link{qri}}, \code{\link{qsr}}, \code{\link{palma_ratio}},
#'   \code{\link{ratio_quantiles}}, \code{\link{rescaled_bootstrap}}
#'
#' @family inequality indicators based on quantiles
#'
#' @importFrom stats setNames
#'
#' @examples
#' data(synthouse)
#' eq <- synthouse$eq_income
#' w  <- synthouse$weight
#'
#' # Point estimates only
#' inequantiles(eq, weights = w)
#'
#' # Subset of indicators
#' inequantiles(eq, weights = w, indicators = c("qri", "palma"))
#'
#' # Custom percentile ratio (P80/P20 instead of P90/P10)
#' inequantiles(eq, weights = w, indicators = "ratio_quantiles",
#'              prob_num = 0.80, prob_den = 0.20)
#'
#' \donttest{
#' # With bootstrap standard errors (complex design)
#' inequantiles(eq, weights = w,
#'              se = TRUE, B = 100, seed = 42,
#'              data = synthouse, strata = "NUTS2",
#'              psu = "municipality")
#' }
#'
#' @export

inequantiles <- function(y,
                         weights    = NULL,
                         indicators = "all",
                         se         = FALSE,
                         type       = 6,
                         na.rm      = TRUE,
                         M          = 100,
                         prob_num   = 0.90,
                         prob_den   = 0.10,
                         B          = 200,
                         seed       = NULL,
                         data       = NULL,
                         strata     = NULL,
                         psu        = NULL,
                         N_h        = NULL,
                         m_h        = NULL) {

  # =========================================================================
  # VALIDATE indicators
  # =========================================================================

  valid <- c("qri", "qsr", "palma", "ratio_quantiles")
  if (identical(indicators, "all")) indicators <- valid

  unknown <- setdiff(indicators, valid)
  if (length(unknown) > 0)
    stop("Unknown indicator(s): ", paste(unknown, collapse = ", "),
         ". Valid choices: ", paste(valid, collapse = ", "), ".")

  # Label for percentile ratio (e.g. "p90p10")
  rq_name <- paste0("p", round(prob_num * 100), "p", round(prob_den * 100))

  # =========================================================================
  # POINT ESTIMATES
  # =========================================================================

  estimates <- c(
    if ("qri" %in% indicators)
      c(qri = qri(y, weights = weights, M = M, type = type, na.rm = na.rm)),
    if ("qsr" %in% indicators)
      c(qsr = qsr(y, weights = weights, type = type, na.rm = na.rm)),
    if ("palma" %in% indicators)
      c(palma = palma_ratio(y, weights = weights, type = type, na.rm = na.rm)),
    if ("ratio_quantiles" %in% indicators)
      setNames(
        ratio_quantiles(y, weights = weights,
                        prob_numerator   = prob_num,
                        prob_denominator = prob_den,
                        type = type, na.rm = na.rm),
        rq_name
      )
  )

  # =========================================================================
  # STANDARD ERRORS (single bootstrap loop over all indicators)
  # =========================================================================

  se_values <- NULL

  if (se) {
    if (is.null(data))
      stop("'data' must be provided when se = TRUE.")
    if (is.null(strata))
      stop("'strata' must be provided when se = TRUE.")

    # Attach y and weights to a copy of data so rescaled_bootstrap can find them
    data$.y_ineq <- y
    y_col <- ".y_ineq"

    if (!is.null(weights)) {
      data$.w_ineq <- weights
      w_col <- ".w_ineq"
    } else {
      w_col <- NULL
    }

    # Evaluate each indicator on each bootstrap replicate
    multi_est <- function(y_b, w_b = NULL) {
      c(
        if ("qri" %in% indicators)
          c(qri = qri(y_b, weights = w_b, M = M, type = type, na.rm = na.rm)),
        if ("qsr" %in% indicators)
          c(qsr = qsr(y_b, weights = w_b, type = type, na.rm = na.rm)),
        if ("palma" %in% indicators)
          c(palma = palma_ratio(y_b, weights = w_b, type = type, na.rm = na.rm)),
        if ("ratio_quantiles" %in% indicators)
          setNames(
            ratio_quantiles(y_b, weights = w_b,
                            prob_numerator   = prob_num,
                            prob_denominator = prob_den,
                            type = type, na.rm = na.rm),
            rq_name
          )
      )
    }

    boot_result <- rescaled_bootstrap(
      data      = data,
      y         = y_col,
      strata    = strata,
      psu       = psu,
      weights   = w_col,
      N_h       = N_h,
      m_h       = m_h,
      estimator = multi_est,
      by_strata = FALSE,
      B         = B,
      seed      = seed
    )

    se_values <- sqrt(boot_result$variance)
    design <- boot_result$design
  }

  # =========================================================================
  # OUTPUT
  # =========================================================================

  result <- list(
    estimates  = estimates,
    se         = se_values,
    B          = if (se) B else NULL,
    design     = if (se) design else NULL,
    call       = match.call()
  )

  class(result) <- "inequantiles"
  result
}


#' @export
print.inequantiles <- function(x, digits = 4, ...) {
  cat("Quantile-based inequality indicators\n")
  cat("-------------------------------------\n")

  if (!is.null(x$se)) {
    out <- data.frame(
      Estimate = round(x$estimates, digits),
      SE       = round(x$se,        digits),
      row.names = names(x$estimates)
    )
  } else {
    out <- data.frame(
      Estimate  = round(x$estimates, digits),
      row.names = names(x$estimates)
    )
  }

  print(out)

  if (!is.null(x$B))
    cat(sprintf("\nBootstrap replicates: %d\n", x$B))

  invisible(x)
}
