#' Quantile-based inequality indicators
#'
#' Estimates one or more quantile-based inequality indicators simultaneously —
#' QRI, quantile-based share ratio (QSR, Palma, or custom), percentile ratio —
#' together with the Gini coefficient as a widely used benchmark. When standard
#' errors are requested, all indicators are evaluated on the same bootstrap
#' replicates, ensuring full comparability.
#'
#' @param y A numeric vector of strictly positive values (e.g. income, wealth, expenditure).
#' @param weights A numeric vector of sampling weights. If \code{NULL},
#'   all observations are equally weighted.
#' @param indicators Character vector specifying which indicators to compute.
#'   Use \code{"all"} (default) for all, or any subset of
#'   \code{"qri"}, \code{"qsr"}, \code{"palma"}, \code{"ratio_quantiles"}, \code{"gini"}.
#'   \code{"qsr"} (quintile share ratio) and \code{"palma"} (Palma index) are special cases of \code{\link{share_ratio}}.
#'   \code{"ratio_quantiles"} computes the P90/P10.
#' @param se Logical; if \code{TRUE}, standard errors are estimated via the
#'   rescaled bootstrap; see \code{\link{rescaled_bootstrap}}. Requires \code{data} and \code{strata} (see below).
#' @param type Quantile estimation type: integer \code{4}--\code{9} or
#'   \code{"HD"} for Harrell--Davis (default: \code{6}). See \code{\link{csquantile}}.
#' @param na.rm Logical; remove missing values before computing? Default:
#'   \code{TRUE}.
#' @param M Integer; number of quantile-ratio grid points for the QRI
#'   (default: \code{100}). Only used when the QRI is estimated; see \code{\link{qri}}.
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
#' @param verbose Logical; if \code{TRUE} (default), displays a progress bar during bootstrap iterations.
#'
#' @returns A list with components:
#'   \item{estimates}{Numeric vector of point estimates of inequality indicators.}
#'   \item{se}{Numeric vector of standard errors, or \code{NULL} when
#'     \code{se = FALSE}.}
#'   \item{B}{Number of bootstrap replicates used, or \code{NULL}.}
#'   \item{design}{Sampling design type detected by the bootstrap, or \code{NULL} when \code{se = FALSE}.}
#'   \item{call}{The matched function call.}
#'
#' @details
#' All quantile-based indicators are computed from the same specified
#' \code{\link{csquantile}} type. When \code{se = TRUE}, a \emph{single}
#' bootstrap loop is run through the rescaled bootstrap method
#' (see \code{\link{rescaled_bootstrap}}) and all indicators are evaluated on
#' each replicate, so standard errors are based on identical resamples and are
#' directly comparable.
#'
#' The Gini coefficient is estimated following \insertCite{langel2013variance}{inequantiles},
#' equation 6, using a weighted formula based on cumulative weight sums.
#'
#'
#' @seealso \code{\link{qri}}, \code{\link{share_ratio}},
#'   \code{\link{ratio_quantiles}}, \code{\link{rescaled_bootstrap}}
#'
#' @family inequality indicators based on quantiles
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
#'
#' \donttest{
#' # With bootstrap standard errors (complex design)
#' inequantiles(eq, weights = w,
#'              se = TRUE, B = 50, seed = 42,
#'              data = synthouse, strata = "NUTS2",
#'              psu = "municipality",
#'              verbose = FALSE)
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
                         B          = 200,
                         seed       = NULL,
                         data       = NULL,
                         strata     = NULL,
                         psu        = NULL,
                         N_h        = NULL,
                         m_h        = NULL,
                         verbose    = TRUE) {

  # =========================================================================
  # VALIDATE indicators
  # =========================================================================

  valid <- c("qri", "qsr", "palma", "ratio_quantiles", "gini")
  if (identical(indicators, "all")) indicators <- valid

  unknown <- setdiff(indicators, valid)
  if (length(unknown) > 0)
    stop("Unknown indicator(s): ", paste(unknown, collapse = ", "),
         ". Valid choices: ", paste(valid, collapse = ", "), ".")


  # =========================================================================
  # POINT ESTIMATES
  # =========================================================================

  estimates <- c(
    if ("qri" %in% indicators)
      c(qri = qri(y, weights = weights, M = M, type = type, na.rm = na.rm)),
    if ("qsr" %in% indicators)
      c(qsr = share_ratio(y, weights = weights, type = type, na.rm = na.rm,
                          prob_numerator = 0.80, prob_denominator = 0.20)),
    if ("palma" %in% indicators)
      c(palma = share_ratio(y, weights = weights, type = type, na.rm = na.rm,
                            prob_numerator = 0.90, prob_denominator = 0.40)),
    if ("ratio_quantiles" %in% indicators)
     c(ratio_quantiles(y, weights = weights,
                        prob_numerator   = 0.90,
                        prob_denominator = 0.10,
                        type = type, na.rm = na.rm)),
    if ("gini" %in% indicators)
      c(gini = .gini_coef(y, weights))
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
          c(qsr = share_ratio(y_b, weights = w_b, type = type, na.rm = na.rm,
                              prob_numerator = 0.80, prob_denominator = 0.20)),
        if ("palma" %in% indicators)
          c(palma = share_ratio(y_b, weights = w_b, type = type, na.rm = na.rm,
                                prob_numerator = 0.90, prob_denominator = 0.40)),
        if ("ratio_quantiles" %in% indicators)
          c(ratio_quantiles(y_b, weights = w_b,
                            prob_numerator   = 0.90,
                            prob_denominator = 0.10,
                            type = type, na.rm = na.rm)),
        if ("gini" %in% indicators)
          c(gini = .gini_coef(y_b, w_b))
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
      seed      = seed,
      verbose   = verbose
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


#' @param x An object of class \code{"inequantiles"}.
#' @param digits Integer; number of decimal places for rounding (default: \code{4}).
#' @param ... Further arguments passed to or from other methods.
#' @returns The argument \code{x}, invisibly.
#' @rdname inequantiles
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
