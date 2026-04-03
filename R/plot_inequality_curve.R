# Internal environment to accumulate legend entries across overlaid curves.
# Reset each time a fresh plot is started (add = FALSE).
.ineq_legend <- new.env(parent = emptyenv())
.ineq_legend$entries <- list()

#' Plot the Inequality Curve
#'
#' Plots the inequality curve \eqn{R(p) = Q(p/2) / Q(1 - p/2)} over
#' \eqn{p \in [0, 1]}, from either sampling survey data or a parametric
#' distribution. The shaded area between the curve and the line \eqn{R(p) = 1}
#' equals the QRI.
#'
#' @param y Numeric vector of strictly positive values (e.g. income). Provide
#'   either \code{y} (empirical mode) or \code{qfunction} (parametric mode),
#'   not both.
#' @param qfunction A parametric quantile function, e.g. \code{qlnorm},
#'   \code{qweibull}. Provide either \code{qfunction} or \code{y}, not both.
#' @param qfun_args Named list of additional arguments passed to
#'   \code{qfunction} (e.g. if \code{qfunction} == \code{qlnorm}, then provide
#'   \code{list(meanlog = 9, sdlog = 0.6)}).
#' @param weights Numeric vector of sampling weights. Only used in estimation
#'   mode. If \code{NULL}, all observations are equally weighted.
#' @param M Integer; number of grid points for evaluating \eqn{R(p)}
#'   (default: \code{200}).
#' @param type Quantile estimation type (integer 4--9 or \code{"HD"}).
#'   Only used in empirical mode. Default: \code{6}.
#' @param na.rm Logical; remove missing values? Default: \code{TRUE}.
#' @param shade Logical; if \code{TRUE} (default), shades the area between
#'   \eqn{R(p)} and the equality line \eqn{R(p) = 1}. The shaded area equals
#'   the QRI.
#' @param add Logical; if \code{TRUE}, adds the curve to an existing plot
#'   without redrawing axes. Each successive call appends a new entry to the
#'   legend on the right of the plot.
#' @param col Colour of the inequality curve (default: \code{"steelblue"}).
#' @param shade_col Colour for the shaded area. Defaults to a transparent
#'   version of \code{col}.
#' @param lwd Line width (default: \code{1.5}).
#' @param lty Line type (default: \code{1}).
#' @param legend_qri Logical; if \code{TRUE} (default), maintains a legend
#'   outside the plot area on the right, showing the QRI for each curve. When
#'   \code{add = TRUE} the new curve is appended below the previous entries.
#' @param label Character string; overrides the auto-generated legend label
#'   (\code{"QRI = X.XXXX"}). Useful for giving curves descriptive names.
#'   Default: \code{NULL} (auto-label).
#' @param xlab x-axis label (default: \code{"p"}).
#' @param ylab y-axis label (default: \code{"R(p)"}).
#' @param main Plot title (default: \code{"Inequality curve"}).
#'
#' @returns Beyond the plot, a named list with three elements:
#'   \item{p}{Numeric vector of grid points in \eqn{[0, 1]}.}
#'   \item{Rp}{Numeric vector of \eqn{R(p)} values at each grid point.}
#'   \item{qri}{The estimated QRI (area between the equality line and the curve).}
#'   The list is returned invisibly, meaning it is not printed to the console
#'   when the function is called without assignment. Assign the output to a
#'   variable (e.g. \code{out <- plot_inequality_curve(...)}) to inspect it.
#'
#' @details
#' The inequality curve \eqn{R(p)} plots the ratio of symmetric quantiles
#' around the median:
#' \deqn{R(p) = \frac{Q(p/2)}{Q(1 - p/2)}, \quad p \in [0, 1]}, against \eqn{p}.
#' For a perfectly equal distribution \eqn{R(p) = 1} for all \eqn{p}, and the
#' curve coincides with the horizontal line at 1. The further the curve lies
#' below the equality line, the more unequal the distribution. The QRI is the
#' area between the equality line and the curve.
#'
#' Boundary values \eqn{R(0) = 0} and \eqn{R(1) = 1} are set by convention
#' (see \insertCite{prendergast2018simple;textual}{inequantiles}).
#'
#' Multiple curves can be overlaid by calling the function repeatedly with
#' \code{add = TRUE}. The legend outside the plot accumulates an entry for
#' each curve automatically.
#'
#' @references
#' \insertRef{prendergast2018simple}{inequantiles}
#'
#' \insertRef{scarpa2025inference}{inequantiles}
#'
#' @seealso \code{\link{qri}} for the sample-based QRI estimator,
#'   \code{\link{superpop_qri}} for the theoretical QRI of a parametric
#'   distribution.
#'
#' @family inequality indicators based on quantiles
#'
#' @importFrom grDevices adjustcolor
#' @importFrom graphics abline legend lines par plot polygon
#'
#' @examples
#' # -----------------------------------------------------------------
#' # Parametric mode: single curve
#' # -----------------------------------------------------------------
#' plot_inequality_curve(
#'   qfunction = qlnorm,
#'   qfun_args = list(meanlog = 9, sdlog = 0.9),
#'   main = "Log-Normal inequality curve"
#' )
#'
#' # -----------------------------------------------------------------
#' # Overlay multiple curves — legend accumulates automatically
#' # -----------------------------------------------------------------
#' plot_inequality_curve(
#'   qfunction = qlnorm,
#'   qfun_args = list(meanlog = 9, sdlog = 0.3),
#'   main  = "Log-Normal inequality curves",
#'   col   = "steelblue",
#'   label = "LogN(9, 0.3)"
#' )
#' plot_inequality_curve(
#'   qfunction = qlnorm,
#'   qfun_args = list(meanlog = 9, sdlog = 0.9),
#'   col   = "tomato", lty = 2, add = TRUE,
#'   label = "LogN(9, 0.9)"
#' )
#'
#' # -----------------------------------------------------------------
#' # Empirical mode: survey data with sampling weights
#' # -----------------------------------------------------------------
#' data(synthouse)
#' out <- plot_inequality_curve(
#'   y       = synthouse$eq_income,
#'   weights = synthouse$weight,
#'   main    = "Inequality curve — synthouse"
#' )
#'
#' # Inspect the returned list
#' out$qri          # estimated QRI
#' head(out$p)      # grid points
#' head(out$Rp)     # R(p) values
#'
#' @export
plot_inequality_curve <- function(y          = NULL,
                                  qfunction  = NULL,
                                  qfun_args  = list(),
                                  weights    = NULL,
                                  M          = 200,
                                  type       = 6,
                                  na.rm      = TRUE,
                                  shade      = TRUE,
                                  add        = FALSE,
                                  col        = "steelblue",
                                  shade_col  = NULL,
                                  lwd        = 1.5,
                                  lty        = 1,
                                  legend_qri = TRUE,
                                  label      = NULL,
                                  xlab       = "p",
                                  ylab       = "R(p)",
                                  main       = "Inequality curve") {

  # -----------------------------------------------------------------------
  # Validate
  # -----------------------------------------------------------------------
  if (!is.null(y) && !is.null(qfunction))
    stop("Provide either 'y' (empirical) or 'qfunction' (parametric), not both.")
  if (is.null(y) && is.null(qfunction))
    stop("Provide either 'y' (numeric vector) or 'qfunction' (quantile function).")

  # -----------------------------------------------------------------------
  # Grid: boundary points 0 and 1 set by convention (R(0)=0, R(1)=1)
  # -----------------------------------------------------------------------
  p_mid <- ((1:M) - 0.5) / M
  p     <- c(0, p_mid, 1)

  # -----------------------------------------------------------------------
  # Compute R(p)
  # -----------------------------------------------------------------------
  if (!is.null(y)) {
    q_lower <- csquantile(y, weights = weights, probs = p_mid / 2,
                          type = type, na.rm = na.rm)
    q_upper <- csquantile(y, weights = weights, probs = 1 - p_mid / 2,
                          type = type, na.rm = na.rm)
    Rp_mid <- q_lower / q_upper

    est <- mean(1 - Rp_mid, na.rm = TRUE)
  } else {
    Rp_mid <- do.call(qfunction, c(list(p_mid / 2),     qfun_args)) /
              do.call(qfunction, c(list(1 - p_mid / 2), qfun_args))

    est <- do.call(superpop_qri, c(list(qfunction), qfun_args))
  }

  Rp      <- c(0, Rp_mid, 1)
  qri_val <- round(est, digits = 2)

  if (is.null(shade_col))
    shade_col <- adjustcolor(col, alpha.f = 0.2)

  # -----------------------------------------------------------------------
  # Fresh plot: reset legend state and widen right margin for legend
  # -----------------------------------------------------------------------
  if (!add) {
    .ineq_legend$entries <- list()
    if (legend_qri) {
      old_mar      <- par("mar")
      new_mar      <- old_mar
      new_mar[4]   <- max(old_mar[4], 10)
      par(mar = new_mar)
    }
    plot(p, Rp,
         type = "n",
         xlim = c(0, 1), ylim = c(0, 1),
         xlab = xlab, ylab = ylab, main = main,
         xaxs = "i", yaxs = "i", bty = "l")
    abline(h = 1, lty = 2, col = "grey50")
  }

  # -----------------------------------------------------------------------
  # Shade and curve
  # -----------------------------------------------------------------------
  if (shade) {
    valid <- !is.nan(Rp)
    polygon(c(p[valid], rev(p[valid])),
            c(Rp[valid], rep(1, sum(valid))),
            col = shade_col, border = NA)
  }

  lines(p, Rp, col = col, lwd = lwd, lty = lty)

  # -----------------------------------------------------------------------
  # Legend: accumulate entries, redraw full legend outside plot on the right
  # -----------------------------------------------------------------------
  if (legend_qri) {
    leg_text <- if (is.null(label)) paste0("QRI = ", sprintf("%.2f", qri_val)) else label
    .ineq_legend$entries <- c(
      .ineq_legend$entries,
      list(list(text = leg_text, col = col, lwd = lwd, lty = lty))
    )

    all_text <- sapply(.ineq_legend$entries, `[[`, "text")
    all_col  <- sapply(.ineq_legend$entries, `[[`, "col")
    all_lwd  <- sapply(.ineq_legend$entries, `[[`, "lwd")
    all_lty  <- sapply(.ineq_legend$entries, `[[`, "lty")

    op <- par(xpd = NA)
    on.exit(par(op), add = TRUE)
    usr <- par("usr")
    legend(x = usr[2] + 0.03 * diff(usr[1:2]),
           y = usr[4],
           legend = all_text,
           col    = all_col,
           lwd    = all_lwd,
           lty    = all_lty,
           bty    = "n",
           xjust  = 0,
           yjust  = 1)
  }

  invisible(list(p = p, Rp = Rp, qri = qri_val))
}
