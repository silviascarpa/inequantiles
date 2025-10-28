#' Gini Coefficient for Grouped Data
#'
#' Computes the Gini coefficient from grouped income data based on linear
#' interpolation of income shares.
#'
#' @param Y Numeric vector of total amounts per group.
#'   Each element represents the sum of the variable of interest (e.g. income, wealth,
#'   consumption) for all individuals or units in that group.
#' @param freq Numeric vector of frequencies (number of individuals or units) per group.
#'   Must have the same length as \code{Y}.
#'
#' @return A numeric value representing the estimated Gini coefficient on grouped data.
#'   The Gini coefficient ranges from 0 (perfect equality) to 1 (complete inequality).
#'   Note that it assumes equality within groups.
#'
#' @details
#' Consider grouped data divided into \eqn{L} classes with known boundaries,
#' observed frequencies \eqn{f_1, \ldots, f_L} and total amounts \eqn{Y_1, \ldots, Y_J}.
#' The Gini coefficient is approximated by linear interpolation of cumulative shares, as:
#'
#' \deqn{G(s_j, u_j) \approx 1 - \sum_{j=1}^{J} (s_j + s_{j-1})(u_j - u_{j-1})}
#'
#' where:
#' \itemize{
#'   \item \eqn{p_j = f_j / \sum_{i=1}^{J} f_i} is the population share of group \eqn{j};
#'   \item \eqn{c_j = Y_j / \sum_{i=1}^{J} Y_i} is the share of the variable of interest in group \eqn{j};
#'   \item \eqn{s_j = \sum_{k=1}^{j} c_k} is the cumulative share of the variable up to group \eqn{j};
#'   \item \eqn{u_j = \sum_{k=1}^{j} p_k} is the cumulative population share up to group \eqn{j};
#'   \item \eqn{s_0 = u_0 = 0} by convention.
#' }
#'
#' This formula computes twice the area between the egalitarian line (perfect equality)
#' and the Lorenz curve obtained by linearly interpolating the points \eqn{(u_j, s_j)}.
#' Since it assumes all observations within a group have identical values, it provides
#' a \emph{lower-bound} estimate of the true Gini coefficient -- actual inequality may be larger
#' \insertCite{jorda2021inequality}{inequantiles}.
#'
#'
#' @note
#' \strong{Important limitations:}
#' \itemize{
#'   \item This is a \strong{lower bound} approximation. The true Gini coefficient
#'     may be higher due to within-group inequality.
#'   \item The bias magnitude depends on the number of groups and how they are defined.
#'   \item Groups should ideally be defined to maximize within-group homogeneity.
#' }
#'
#'
#' @references
#' \insertRef{jorda2021inequality}{inequantiles}
#'
#' @examples
#' income_freq <- c(1200, 1800, 1500, 800, 400, 20, 10)
#' income_tot <- c(18800, 16300, 44700, 33900, 21500, 22100, 98300)
#'
#' gini_grouped(Y = income_tot, freq = income_freq)
#'
#' @seealso
#' \code{\link{qri_grouped}} for computing the quantile ratio index from grouped data.
#'
#' @export
gini_grouped <- function(Y, freq) {
  # Input validation
  if (length(Y) != length(freq)) {
    stop("Y and freq must have the same length")
  }

  if (any(freq < 0)) {
    stop("freq must be non-negative")
  }

  if (any(Y < 0)) {
    warning("Y must be non-negative")
  }


  if (sum(freq) == 0) {
    stop("Total population (sum of freq) cannot be zero")
  }

  # Ensure numeric input
  Y <- as.numeric(Y)  # Total per group
  freq <- as.numeric(freq)  # Number of people per group

  # Total population and total variable
  Tot.ppl <- sum(freq)
  Tot.Y <- sum(Y)

  if (Tot.Y == 0) {
    warning("Total income is zero. Returning Gini = 0")
    return(0)
  }

  # Compute income and population shares
  c_j <- Y / Tot.Y  # Income share per group
  p_j <- freq / Tot.ppl  # Population share per group

  # Compute cumulative proportions
  s_j <- cumsum(c_j)   # Cumulative income share
  u_j <- cumsum(p_j)   # Cumulative population share

  # Lagged cumulative shares (starting from 0)
  s_jLag <- c(0, s_j[-length(s_j)])
  u_jLag <- c(0, u_j[-length(u_j)])

  # Compute Gini coefficient using equation (1) from Jordá et al. (2020)
  # G ≈ 1 - Σ(s_j + s_{j-1})(u_j - u_{j-1})
  G <- 1 - sum((s_j + s_jLag) * (u_j - u_jLag))

  return(G)
}
