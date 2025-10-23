#' Rescaled Bootstrap Variance Estimation
#'
#'
#' Implements the rescaled bootstrap method for variance estimation in survey data,
#' supporting both stratified simple random sampling and two-stage complex designs.
#' Based on Rao & Wu (1988, 1992).
#'
#' @param data A data frame containing the survey data
#' @param estimator A function that takes the data (and weights if needed)
#'   and returns the estimate. For simple designs: function(data).
#'   For complex designs: function(data, weights)
#' @param strata A character string or vector specifying the stratification variable
#' @param psu A character string specifying the Primary Sampling Units variable (required for a multistage sampling design)
#' @param weights A character string specifying the weight variable (required if design = "complex")
#' @param N_h Optional vector of stratum population size, needed for computing the  finite population correction (for stratified srs).
#'            Can be a vector of length H (one per stratum) or a single value for all strata.
#'           Default is NULL (no FPC applied)
#' @param B Integer, number of bootstrap replicates (default = 200)
#' @param m_h Optional vector of bootstrap sample sizes per stratum. If NULL, computed as floor((n_h-2)^2/(n_h-1))
#' @param seed Optional integer for reproducibility
#'
#'
#' @return A list containing:
#'   \item{variance}{Bootstrap variance estimate}
#'   \item{boot_estimates}{Vector of B bootstrap estimates}
#'   \item{B}{Number of bootstrap replicates}
#'   \item{by_strata}{if variance is computed by stratum or not}
#'   \item{design}{if the sampling design is complex (with sampling weights) or simple}
#'   \item{strata_info}{Returns information about number of observations/PSUs per stratum}
#'   \item{call}
#'
#'
#'#'#' @references
#' Rao, J.N.K. & Wu, C.F.J. (1988). Resampling inference with complex survey data.
#'   JASA, 83(401), 231-241.
#' Rao, J.N.K., Wu, C.F.J. & Yue, K. (1992). Some recent work on resampling
#'   methods for complex surveys. Survey Methodology, 18(2), 209-217.
#'  Kolenikov, S. (2010). Resampling variance estimation for complex survey data.
#'    The Stata Journal, 10(2), 165-199.
#'  Scarpa, Silvia, Maria Rosaria Ferrante, and Stefan Sperlich.
#'    "INFERENCE FOR THE QUANTILE RATIO INEQUALITY INDEX IN THE CONTEXT OF SURVEY DATA.",
#'    Journal of Survey Statistics and Methodology (2025): smaf024.
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




