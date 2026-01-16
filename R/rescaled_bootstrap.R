#' Rescaled Bootstrap Variance Estimation
#'
#'
#' Implements the *rescaled bootstrap* method for variance estimation in survey data,
#' supporting both stratified simple random sampling and multistage complex designs.
#'
#' @param data A data frame containing the survey data.
#' @param y A character string specifying the variable name to be used for the target variable.
#' @param strata A character string specifying the stratification variable.
#' @param N_h Optional vector of stratum population sizes, used for the finite population correction (FPC).
#'   Can be a single value (applied to all strata) or one value per stratum.
#' @param psu Optional character string specifying the Primary Sampling Unit (PSU) variable.
#'   Required for multistage complex designs.
#' @param weights Optional character string specifying the sampling weight variable.
#'   Required for complex designs with unequal inclusion probabilities.
#' @param estimator A function that computes the statistic of interest, accepting arguments
#'   \code{estimator(y, weights)} for complex designs or \code{estimator(y)} for simple designs.
#' @param by_strata Logical; if \code{TRUE}, variances are computed separately by stratum.
#' @param B Integer; number of bootstrap replicates (default = 200).
#' @param m_h Optional vector of bootstrap sample sizes per stratum (PSUs for complex designs).
#'   If \code{NULL}, defaults to \eqn{m_h = \lfloor (n_h - 2)^2 / (n_h - 1) \rfloor}.
#' @param seed Optional integer for reproducibility.
#'
#'
#' @return A list containing:
#'   \item{variance}{Bootstrap variance estimate}
#'   \item{boot_estimates}{Vector of B bootstrap estimates}
#'   \item{B}{Number of bootstrap replicates}
#'   \item{by_strata}{if variance is computed by stratum or overall}
#'   \item{design}{if the sampling design is complex (with sampling weights) or simple}
#'   \item{strata_info}{Returns information about number of observations/PSUs per stratum}
#'   \item{call}{The matched function call.}
#'
#' @details
#' The rescaled bootstrap is a resampling technique designed for complex survey data that preserves
#' stratification and primary sampling unit (PSU) structure, providing consistent variance estimation
#' for both smooth and non-smooth statistics.
#' The methodology is based on \insertCite{rao1988resampling;textual}{inequantiles} and
#' \insertCite{rao1992some;textual}{inequantiles}.
#'
#'
#' \strong{(1) Stratified Simple Random Sampling}
#'
#' Consider a finite population divided into \eqn{H} strata, each of size \eqn{N_h}, with a sample of size \eqn{n_h}
#' selected independently in each stratum. For each \eqn{b} bootstrap replicate, \eqn{b = \ldots, B}:
#' \enumerate{
#'   \item Draw a bootstrap sample of size \eqn{m_h} with replacement from the \eqn{n_h} sampled units.
#'         By default, \eqn{m_h = \lfloor (n_h - 2)^2 / (n_h - 1) \rfloor \approx n_h - 3}.
#'   \item Compute rescaled bootstrap values:
#'         \deqn{
#'           \tilde{y}_{hj}^{*(b)} = \bar{y}_h +
#'           \sqrt{\frac{m_h(1-f_h)}{n_h - 1}}
#'           (y_{hj}^{*(b)} - \bar{y}_h), }
#'         where \eqn{y_{hj}^*} is the bootstrap observation \eqn{f_h = n_h / N_h} is the sampling fraction and \eqn{\bar{y}_h} is the sample stratum mean.
#'   \item Compute the statistic of interest \eqn{\hat{\theta}^{*(b)}} using rescaled values.
#' }
#'
#' The bootstrap variance is then given by:
#' \deqn{
#'   \widehat{V}_{boot}(\hat{\theta}) =
#'   \frac{1}{B-1} \sum_{b=1}^{B}
#'   \left( \hat{\theta}^{*(b)} - \bar{\theta}^{*} \right)^2,
#'   \qquad
#'   \bar{\theta}^{*} = \frac{1}{B} \sum_{b=1}^{B} \hat{\theta}^{*(b)}.
#' }
#'
#' \strong{(2) Two-Stage Stratified Sampling (Complex Designs)}
#'
#' For designs with PSUs and sampling weights:
#' \enumerate{
#'   \item Within each stratum \eqn{h}, draw \eqn{m_h} PSUs with replacement from the \eqn{n_h} sampled PSUs.
#'         By default, \eqn{m_h = \lfloor (n_h - 2)^2 / (n_h - 1) \rfloor \approx n_h - 3}.
#'   \item Let \eqn{m_{hi}^{(b)}} denote the number of times PSU \eqn{i} is selected in replicate \eqn{b}.
#'         Each observation in the \eqn{i}-th PSU is assigned a rescaled bootstrap weight:
#'         \deqn{
#'           w_{hij}^{*(b)} =
#'           \left[
#'             1 - c_h + c_h \frac{n_h}{m_h} m_{hi}^{(b)}
#'           \right] w_{hij},
#'           \qquad
#'           c_h = \sqrt{\frac{m_h}{n_h - 1}}.
#'         }
#'         \eqn{w_{hij}} is the sampling weight associtaed to individual
#'             \eqn{j} in PSU \eqn{i} in stratum \eqn{h}-
#'   \item The statistic \eqn{\hat{\theta}^{*(b)}} is computed using the rescaled weights.
#' }
#'
#' The rescaled bootstrap variance estimate is then:
#' \deqn{
#'   \widehat{V}_{boot}(\hat{\theta}) =
#'   \frac{1}{B-1} \sum_{b=1}^{B}
#'   \left( \hat{\theta}^{*(b)} - \bar{\theta}^{*} \right)^2.
#' }
#'
#' @examples
#' \donttest{
#' data(synthouse)
#'
#'
#' # ================================================================
#' # Example 1: Stratified Simple Random Sampling (SRS)
#' # ================================================================
#'
#'
#' # Use NUTS2 as strata
#' set.seed(123)
#'
#' # Simulate population sizes per stratum (for FPC)
#' N_values <- sample(2000:5000, length(unique(synthouse$NUTS2)), replace = TRUE)
#' names(N_values) <- sort(unique(synthouse$NUTS2))
#'
#' # Define a simple mean estimator
#' mean_estimator <- function(y) mean(y, na.rm = TRUE)
#'
#' # Apply the rescaled bootstrap under stratified SRS
#' boot_srs <- rescaled_bootstrap(
#'   data = synthouse,
#'   y = "eq_income",
#'   strata = "NUTS2",
#'   N_h = N_values,
#'   estimator = mean_estimator,
#'   by_strata = TRUE,
#'   B = 50,  # small number for illustration
#'   seed = 123
#' )
#'
#' # View results
#' boot_srs$variance
#'
#'
#' # ================================================================
#' # Example 2: Two-stage Complex Design
#' # ================================================================
#'
#' # PSU = municipality, Strata = NUTS2, weights = weight, y = eq_income
#'
#' # Define a weighted mean estimator
#' wmean_estimator <- function(y, weights) {
#'   sum(y * weights, na.rm = TRUE) / sum(weights, na.rm = TRUE)
#' }
#'
#' boot_complex <- rescaled_bootstrap(
#'   data = synthouse,
#'   y = "eq_income",
#'   strata = "NUTS2",
#'   psu = "municipality",
#'   weights = "weight",
#'   estimator = wmean_estimator,
#'   by_strata = TRUE,
#'   B = 50,
#'   seed = 456
#' )
#'
#' # Display variance and bootstrap estimates
#' summary(boot_complex$variance)
#'
#' # Strata and PSU summary
#' boot_complex$strata_info
#' }
#'
#'
#' # ================================================================
#' # Note:
#' # These examples use small B for speed. For actual analysis,
#' # use B >= 200 for stable estimates.
#' # ================================================================
#
#'
#'
#' @references
#'
#' \insertRef{rao1988resampling}{inequantiles}
#'
#' \insertRef{rao1992some}{inequantiles}
#'
#' \insertRef{kolenikov2010resampling}{inequantiles}
#'
#' \insertRef{scarpa2025inference}{inequantiles}
#'
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#'
#' @export

rescaled_bootstrap <- function(data,
                               y,
                               strata,
                               N_h = NULL,
                               psu = NULL,
                               weights = NULL,
                               estimator,
                               by_strata = TRUE,
                               B = 200,
                               m_h = NULL,
                               seed = NULL) {

  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Detect design type
  is_complex <- !is.null(weights)
  has_psu <- !is.null(psu)

  # Check for problematic configurations
  if (is_complex && !has_psu) {
    stop(paste0(
      "ERROR: 'weights' specified but 'psu' is NULL.\n",
      "For complex design with weights, you MUST specify PSU variable.\n",
      "Either:\n",
      "  1. Add psu = 'your_psu_column' for complex two-stage design, OR\n",
      "  2. Set weights = NULL for simple stratified design"
    ))
  }

  if (!is_complex && has_psu) {
    stop(paste0(
      "WARNING: 'psu' specified but 'weights' is NULL.\n",
      "PSU variable will be IGNORED in simple design.\n",
      "This will produce variance = 0 (no between-PSU variability captured).\n",
      "Did you mean to specify weights = 'your_weight_column' for complex design?"
    ))
  }

  if (!is_complex && !has_psu && is.null(N_h)) {
    stop(paste0(
      "WARNING: Simple design without finite population correction (N_h = NULL).\n",
      "The FPC will be set to 0, resulting in variance = 0.\n",
      "To get meaningful variance estimates, you should specify N_h (population sizes).\n",
      "Example: N_h = c(1000, 1500, 800) for 3 strata"
    ))
  }



  # ========================================================================
  # SETUP: Create standard column names
  # ========================================================================

  data$stratum_var <- as.character(data[[strata]])
  data$y <- as.numeric(data[[y]])

  if (is_complex) {
    # Complex design: with PSUs and weights
    if (is.null(psu)) {
      stop("For complex design (weights provided), 'psu' must be specified")
    }
    data$psu_var <- as.character(data[[psu]])
    data$weight_var <- as.numeric(data[[weights]])
  }

  # ========================================================================
  # EXTRACT STRATA INFORMATION
  # ========================================================================

  strata_levels <- sort(unique(data$stratum_var))
  n_strata <- length(strata_levels)

  if (is_complex) {
    # Complex: n_h = number of PSUs per stratum
    psu_strata <- unique(data[, c("stratum_var", "psu_var")])
    n_h <- numeric(n_strata)
    names(n_h) <- strata_levels

    for (h in strata_levels) {
      n_h[h] <- sum(psu_strata$stratum_var == h)
    }

  } else {
    # Simple: n_h = number of observations per stratum
    n_h <- numeric(n_strata)
    names(n_h) <- strata_levels

    for (h in strata_levels) {
      n_h[h] <- sum(data$stratum_var == h)
    }
  }

  # ========================================================================
  # PROCESS N_h (population sizes for FPC)
  # ========================================================================

  if (!is_complex && !is.null(N_h)) {
    # Validate N_h
    if (length(N_h) == 1) {
      N_h <- rep(N_h, n_strata)
      names(N_h) <- strata_levels
    } else if (length(N_h) != n_strata) {
      stop(sprintf("'N_h' must be length 1 or %d (number of strata)", n_strata))
    }

    # Check that N_h >= n_h
    if (any(N_h < n_h)) {
      problematic_strata <- strata_levels[N_h < n_h]
      stop(sprintf(
        "ERROR: N_h (population size) must be >= n_h (sample size) for all strata.\n",
        "Problematic strata: %s\n",
        "N_h values: %s\n",
        "n_h values: %s",
        paste(problematic_strata, collapse = ", "),
        paste(N_h[problematic_strata], collapse = ", "),
        paste(n_h[problematic_strata], collapse = ", ")
      ))
    }
  }


  # ========================================================================
  # CALCULATE m_h (bootstrap sample sizes)
  # ========================================================================

  if (is.null(m_h)) {
    # Complex: Rao & Wu formula for PSU sampling
    m_h <- ifelse(n_h <= 3,
                  n_h,
                  floor((n_h - 2)^2 / (n_h - 1)))
    m_h <- pmax(m_h, 2)
  } else {
    if (length(m_h) == 1) {
      m_h <- rep(m_h, n_strata)
      names(m_h) <- strata_levels
    } else if (length(m_h) != n_strata) {
      stop(sprintf("'m_h' must be length 1 or %d (number of strata)", n_strata))
    }
  }


  # ========================================================================
  # PREPARE STORAGE
  # ========================================================================

  if (by_strata) {
    boot_estimates <- matrix(NA, nrow = B, ncol = n_strata)
    colnames(boot_estimates) <- strata_levels
  } else {
    boot_estimates <- numeric(B)
  }

  # Progress bar
  pb <- txtProgressBar(min = 0, max = B, style = 3)

  # ========================================================================
  # BOOTSTRAP LOOP
  # ========================================================================

  if (is_complex) {
    # =====================================================================
    # COMPLEX DESIGN: Sample PSUs and rescale WEIGHTS
    # =====================================================================

    for (b in 1:B) {

      # Create frequency table for this bootstrap replicate
      psu_freq <- data.frame(
        stratum_var = character(),
        psu_var = character(),
        freq = integer(),
        stringsAsFactors = FALSE
      )

      l <- 0
      for (h in strata_levels) {
        psu_h <- psu_strata$psu_var[psu_strata$stratum_var == h]
        n_h_val <- n_h[h]
        l <- l + 1
        m_h_val <- m_h[l]

        # Sample m_h PSUs with replacement
        sampled_psu <- sample(psu_h, size = m_h_val, replace = TRUE)

        # Count frequencies
        freq_table <- table(factor(sampled_psu, levels = psu_h))

        # Add to frequency dataframe
        for (psu_id in psu_h) {
          psu_freq <- rbind(psu_freq, data.frame(
            stratum_var = h,
            psu_var = psu_id,
            freq = as.integer(freq_table[psu_id]),
            stringsAsFactors = FALSE
          ))
        }
      }

      # Merge frequencies with original data
      data_boot <- merge(data, psu_freq, by = c("stratum_var", "psu_var"), all.x = TRUE)
      data_boot$freq[is.na(data_boot$freq)] <- 0

      # RESCALE WEIGHTS - Rao & Wu (1992) formula
      data_boot$weight_rescaled <- NA

      l <- 0
      for (h in strata_levels) {
        idx_h <- data_boot$stratum_var == h
        n_h_val <- n_h[h]
        l <- l + 1
        m_h_val <- m_h[l]
        freq_h <- data_boot$freq[idx_h]
        weight_h <- data_boot$weight_var[idx_h]

        c_h <- sqrt(m_h_val / (n_h_val - 1))
        weight_rescaled_h <- ((1 - c_h) + c_h * (n_h_val / m_h_val) * freq_h) * weight_h

        data_boot$weight_rescaled[idx_h] <- weight_rescaled_h
      }

      # Compute bootstrap estimate
      if (by_strata) {
        for (h in strata_levels) {
          data_h <- data_boot[data_boot$stratum_var == h, ]
          est_result <- estimator(data_h$y, data_h$weight_rescaled)
          boot_estimates[b, h] <- est_result
        }
      } else {
        est_result <- estimator(data_boot$y, data_boot$weight_rescaled)
        boot_estimates[b] <- est_result
      }

      setTxtProgressBar(pb, b)
    }

  } else {
    # =====================================================================
    # SIMPLE DESIGN: Sample observations and rescale Y VALUES
    # =====================================================================

    for (b in 1:B) {

      # Container for bootstrap sample
      boot_data_list <- list()
      l = 0
      for (h in strata_levels) {
        data_h <- data[data$stratum_var == h, ]
        n_h_val <- n_h[h]
        l <- l + 1
        m_h_val <- m_h[l]
        N_h_val <- ifelse(is.null(N_h), n_h_val, N_h[l])

        # Sample m_h observations with replacement
        boot_indices <- sample(1:n_h_val, size = m_h_val, replace = TRUE)
        boot_sample_h <- data_h[boot_indices, ]

        # RESCALE Y VALUES - Rao & Wu (1988) formula
        y_bar_h <- mean(data_h$y)
        fpc_h <- 1 - n_h_val / N_h_val # Finite population correction (for SRSWR)

        rescale_factor <- sqrt((m_h_val * fpc_h) / (n_h_val - 1))
        boot_sample_h$y_rescaled <- y_bar_h + rescale_factor * (boot_sample_h$y - y_bar_h)

        boot_data_list[[h]] <- boot_sample_h
      }

      # Combine all strata
      boot_data <- do.call(rbind, boot_data_list)

      # Compute bootstrap estimate
      if (by_strata) {
        for (h in strata_levels) {
          boot_data_h <- boot_data[boot_data$stratum_var == h, ]
          est_result <- estimator(boot_data_h$y_rescaled)
          boot_estimates[b, h] <- est_result
        }
      } else {
        est_result <- estimator(boot_data$y_rescaled)
        boot_estimates[b] <- est_result
      }

      setTxtProgressBar(pb, b)
    }
  }

  close(pb)

  # ========================================================================
  # CALCULATE VARIANCE
  # ========================================================================

  if (by_strata) {
    variance <- apply(boot_estimates, 2, var)
    se <- sqrt(variance)
  } else {
    variance <- var(boot_estimates)
    se <- sqrt(variance)
  }

  # ========================================================================
  # PREPARE OUTPUT
  # ========================================================================

  if (is_complex) {
    strata_info <- data.frame(
      stratum = strata_levels,
      n_psu = n_h,
      m_psu = m_h,
      stringsAsFactors = FALSE
    )
  } else {
    strata_info <- data.frame(
      stratum = strata_levels,
      n_obs = n_h,
      m_obs = m_h,
      stringsAsFactors = FALSE
    )
  }

  result <- list(
    variance = variance,
    boot_estimates = boot_estimates,
    B = B,
    by_strata = by_strata,
    design = ifelse(is_complex, "complex", "simple"),
    strata_info = strata_info,
    call = match.call()
  )

  class(result) <- c("rescaled_bootstrap", "list")

  return(result)
}




