# ==============================================================================
# Generate Synthetic IT-SILC Dataset Based on Empirical Structure
# ==============================================================================
#
# This script generates a realistic synthetic dataset that replicates the
# empirical structure of IT-SILC 2024, using the statistics extracted by
# analyze_itsilc.R
#
# Output: synthouse dataset ready for the package
# ==============================================================================

library(dplyr)
library(tidyr)

set.seed(2024)

# ==============================================================================
# Load IT-SILC structure (from previous analysis)
# ==============================================================================

cat("Loading IT-SILC empirical structure...\n")
itsilc_structure <- readRDS("C:/GitHub/inequantiles/itsilc_structure.rds")

# Target sample size
n_individuals_target <- 20000
n_households_target <- round(n_individuals_target / itsilc_structure$mean_hh_size)

cat("Target sample:\n")
cat("- Households:", n_households_target, "\n")
cat("- Individuals:", n_individuals_target, "\n")
cat("- Avg HH size:", itsilc_structure$mean_hh_size, "\n\n")

# ==============================================================================
# 1. GENERATE HOUSEHOLDS - Sample archetypes
# ==============================================================================

cat("Generating households...\n")

# Sample household archetypes
sampled_archetypes <- sample(
  names(itsilc_structure$archetype_probs),
  size = n_households_target,
  replace = TRUE,
  prob = itsilc_structure$archetype_probs
)

# Create household data frame
households <- data.frame(
  hh_id = sprintf("HH%06d", 1:n_households_target),
  archetype = sampled_archetypes,
  stringsAsFactors = FALSE
)

# Helper function to sample household size
sample_hh_size_default <- function(n_samples) {
  # Get available sizes and their probabilities
  sizes <- as.integer(names(itsilc_structure$hh_size_probs))
  probs <- as.numeric(itsilc_structure$hh_size_probs)

  # If we need sizes not in the distribution, pad with small probabilities
  all_sizes <- 1:5
  if (!all(all_sizes %in% sizes)) {
    missing_sizes <- setdiff(all_sizes, sizes)
    sizes <- c(sizes, missing_sizes)
    probs <- c(probs, rep(0.01, length(missing_sizes)))
    probs <- probs / sum(probs)  # Renormalize
  }

  sample(sizes, n_samples, replace = TRUE, prob = probs)
}

# Determine household characteristics from archetype
households <- households %>%
  mutate(
    hh_size = case_when(
      grepl("Single", archetype) ~ 1L,
      grepl("Couple|Parent_child", archetype) ~ 2L,
      grepl("Family_1child", archetype) ~ 3L,
      grepl("Family_2children", archetype) ~ 4L,
      grepl("Family_3plus", archetype) ~ as.integer(sample(5:7, n(), replace = TRUE)),
      grepl("Multi_adult", archetype) ~ as.integer(sample(3:4, n(), replace = TRUE)),
      TRUE ~ sample_hh_size_default(n())
    ),
    n_children = case_when(
      grepl("1child", archetype) ~ 1L,
      grepl("2children", archetype) ~ 2L,
      grepl("3plus", archetype) ~ as.integer(sample(3:4, n(), replace = TRUE)),
      grepl("Parent_child", archetype) ~ 1L,
      TRUE ~ 0L
    ),
    n_adults = hh_size - n_children
  )

cat("Household archetypes sampled.\n")
cat("Average household size:", mean(households$hh_size), "\n\n")

# ==============================================================================
# 2. GENERATE INDIVIDUALS - Expand households
# ==============================================================================

cat("Generating individuals...\n")

# Expand to individual level
individuals <- households %>%
  slice(rep(1:n(), hh_size)) %>%
  group_by(hh_id) %>%
  mutate(
    position_in_hh = row_number(),
    person_id = sprintf("P%08d", row_number() + (cur_group_id() - 1) * max(hh_size))
  ) %>%
  ungroup()

cat("Total individuals:", nrow(individuals), "\n\n")

# ==============================================================================
# 3. ASSIGN AGE based on position and archetype
# ==============================================================================

cat("Assigning age...\n")

# Helper function to sample age from age_class
sample_age_from_class <- function(age_classes) {
  sapply(age_classes, function(ac) {
    if (is.na(ac)) return(NA)
    age_ranges <- list(
      "0-14" = 0:14,
      "15-17" = 15:17,
      "18-24" = 18:24,
      "25-34" = 25:34,
      "35-49" = 35:49,
      "50-64" = 50:64,
      "65+" = 65:80
    )
    sample(age_ranges[[ac]], 1)
  })
}

# Step 1: Assign age class to head (position 1)
individuals <- individuals %>%
  group_by(hh_id) %>%
  mutate(
    age_class = case_when(
      # HEAD OF HOUSEHOLD (position 1)
      position_in_hh == 1 & grepl("young", archetype) ~ "25-34",
      position_in_hh == 1 & grepl("adult", archetype) ~ "35-49",
      position_in_hh == 1 & grepl("mature", archetype) ~ "50-64",
      position_in_hh == 1 & grepl("elderly", archetype) ~ "65+",
      position_in_hh == 1 & n_children > 0 ~ "35-49",  # Parents with children
      position_in_hh == 1 ~ "35-49",  # Default
      TRUE ~ NA_character_
    )
  ) %>%
  ungroup()

# Generate actual age for heads
individuals <- individuals %>%
  mutate(
    age = ifelse(position_in_hh == 1, sample_age_from_class(age_class), NA_integer_)
  )

# Step 2: Assign age to other household members based on head's age
individuals <- individuals %>%
  group_by(hh_id) %>%
  mutate(
    age_head = first(age),

    age = case_when(
      # Already assigned (head)
      position_in_hh == 1 ~ age,

      # SPOUSE/PARTNER (position 2, adult in couples)
      position_in_hh == 2 & n_adults >= 2 & n_children == 0 & age_head < 35 ~
        as.integer(age_head + sample(-5:5, 1)),  # Similar age
      position_in_hh == 2 & n_adults >= 2 & n_children == 0 & age_head >= 35 & age_head < 65 ~
        as.integer(age_head + sample(-7:7, 1)),  # Similar age
      position_in_hh == 2 & n_adults >= 2 & n_children == 0 & age_head >= 65 ~
        as.integer(age_head + sample(-5:5, 1)),  # Elderly couple

      # SPOUSE/PARTNER in families with children
      position_in_hh == 2 & n_adults >= 2 & n_children > 0 ~
        as.integer(pmax(20, age_head + sample(-5:5, 1))),  # Similar age, min 20

      # CHILDREN (position > n_adults)
      position_in_hh > n_adults & position_in_hh == (n_adults + 1) ~
        as.integer(pmax(0, age_head - sample(20:35, 1))),  # First child: parent 20-35 years older
      position_in_hh > n_adults & position_in_hh == (n_adults + 2) ~
        as.integer(pmax(0, age_head - sample(22:37, 1))),  # Second child
      position_in_hh > n_adults & position_in_hh == (n_adults + 3) ~
        as.integer(pmax(0, age_head - sample(24:39, 1))),  # Third child
      position_in_hh > n_adults ~
        as.integer(pmax(0, age_head - sample(20:40, 1))),  # Other children

      # YOUNG ADULT CHILD (Parent-child household)
      position_in_hh == 2 & grepl("Parent_child", archetype) ~
        as.integer(pmax(18, age_head - sample(20:30, 1))),  # Adult child

      # OTHER ADULTS (multi-adult households)
      position_in_hh <= n_adults ~
        as.integer(pmax(18, age_head + sample(-10:10, 1))),  # Other adults

      TRUE ~ as.integer(sample(20:70, 1))  # Fallback
    )
  ) %>%
  ungroup()

# Step 3: Ensure consistency - children < 18, adjust if needed
individuals <- individuals %>%
  group_by(hh_id) %>%
  mutate(
    # If marked as child but age >= 18, reduce age
    age = ifelse(position_in_hh > n_adults & age >= 18,
                 sample(0:17, 1), age),

    # If parent too young to have children, adjust
    age = ifelse(position_in_hh == 1 & n_children > 0 & age < 20,
                 sample(25:35, 1), age),

    # Ensure age >= 0
    age = pmax(0, age),

    # Recalculate age_class based on actual age
    age_class = cut(age,
                    breaks = c(-1, 14, 17, 24, 34, 49, 64, 100),
                    labels = c("0-14", "15-17", "18-24", "25-34", "35-49", "50-64", "65+"))
  ) %>%
  ungroup()

cat("Age assigned with realistic family structure.\n")
cat("Age distribution:\n")
print(table(individuals$age_class))
cat("\n")

# Check parent-child age differences
parent_child_check <- individuals %>%
  group_by(hh_id) %>%
  filter(n_children > 0) %>%
  summarise(
    age_parent = first(age),
    age_youngest_child = min(age[position_in_hh > n_adults]),
    age_diff = age_parent - age_youngest_child,
    .groups = "drop"
  )

if (nrow(parent_child_check) > 0) {
  cat("Parent-child age difference check:\n")
  cat("Min:", min(parent_child_check$age_diff), "\n")
  cat("Median:", median(parent_child_check$age_diff), "\n")
  cat("Max:", max(parent_child_check$age_diff), "\n\n")
}

# ==============================================================================
# 4. ASSIGN GENDER
# ==============================================================================

cat("Assigning gender...\n")

# First pass: assign gender to all individuals
individuals <- individuals %>%
  mutate(
    gender = sample(1:2, n(), replace = TRUE,
                    prob = itsilc_structure$gender_probs)
  )

# Second pass: for couples, ensure different genders when possible
individuals <- individuals %>%
  group_by(hh_id) %>%
  mutate(
    gender = case_when(
      # For position 2 in couples, try opposite gender of position 1
      position_in_hh == 2 & hh_size >= 2 & n_adults >= 2 ~
        3 - first(gender),  # Opposite of head
      # Keep original for others
      TRUE ~ gender
    )
  ) %>%
  ungroup()

cat("Gender assigned.\n\n")

# ==============================================================================
# 5. ASSIGN EDUCATION based on age
# ==============================================================================

cat("Assigning education level...\n")

# Extract education probabilities by age from structure
edu_by_age_wide <- itsilc_structure$education_by_age %>%
  select(age_class, education_macro, prop) %>%
  pivot_wider(names_from = education_macro, values_from = prop, values_fill = 0)

assign_education <- function(age_class) {
  if (is.na(age_class) || age_class %in% c("0-14", "15-17")) {
    return(NA_character_)
  }

  # Get probabilities for this age class
  probs <- edu_by_age_wide %>% filter(age_class == !!age_class)

  if (nrow(probs) == 0) {
    return(sample(c("Low", "Medium", "High"), 1, prob = c(0.3, 0.5, 0.2)))
  }

  # Extract probabilities, handling missing columns
  prob_high <- if ("High" %in% names(probs)) probs$High else 0
  prob_low <- if ("Low" %in% names(probs)) probs$Low else 0
  prob_medium <- if ("Medium" %in% names(probs)) probs$Medium else 0

  # Ensure we have valid probabilities
  total_prob <- prob_high + prob_low + prob_medium
  if (total_prob == 0) {
    return(sample(c("Low", "Medium", "High"), 1, prob = c(0.3, 0.5, 0.2)))
  }

  # Normalize
  prob_vec <- c(High = prob_high, Low = prob_low, Medium = prob_medium) / total_prob

  sample(names(prob_vec), 1, prob = prob_vec)
}

individuals <- individuals %>%
  rowwise() %>%
  mutate(education_level = assign_education(age_class)) %>%
  ungroup()

cat("Education level assigned.\n")
cat("Education distribution (adults):\n")
print(table(individuals$education_level[!is.na(individuals$education_level)]))
cat("\n")

# ==============================================================================
# 6. ASSIGN OCCUPATION based on age and education
# ==============================================================================

cat("Assigning occupation...\n")

# Extract occupation probabilities
occ_by_age_wide <- itsilc_structure$occupation_by_age %>%
  select(age_class, occupation_macro, prop) %>%
  pivot_wider(names_from = occupation_macro, values_from = prop, values_fill = 0)

assign_occupation <- function(age_class, education) {
  if (is.na(age_class)) {
    return(NA_character_)
  }

  # Get probabilities for this age class
  probs <- occ_by_age_wide %>% filter(age_class == !!age_class)

  if (nrow(probs) == 0) {
    return(sample(c("Employed", "Unemployed", "Other", "Retired", "Student"),
                  1, prob = c(0.5, 0.1, 0.1, 0.2, 0.1)))
  }

  # Extract occupation names and probabilities (excluding age_class column)
  occ_names <- setdiff(names(probs), "age_class")
  occ_probs <- as.numeric(probs[1, occ_names])
  names(occ_probs) <- occ_names

  # Remove zero probabilities to avoid issues
  occ_probs <- occ_probs[occ_probs > 0]

  if (length(occ_probs) == 0) {
    return(sample(c("Employed", "Unemployed", "Other", "Retired", "Student"),
                  1, prob = c(0.5, 0.1, 0.1, 0.2, 0.1)))
  }

  # Adjust for education if available (high education -> more employed)
  if (!is.na(education) && education == "High" && "Employed" %in% names(occ_probs)) {
    occ_probs["Employed"] <- occ_probs["Employed"] * 1.2
    occ_probs <- occ_probs / sum(occ_probs)  # Renormalize
  }

  sample(names(occ_probs), 1, prob = occ_probs)
}

individuals <- individuals %>%
  rowwise() %>%
  mutate(employment_status = assign_occupation(age_class, education_level)) %>%
  ungroup()

cat("Occupation assigned.\n")
cat("Occupation distribution:\n")
print(table(individuals$employment_status))
cat("\n")

# ==============================================================================
# 7. ASSIGN GEOGRAPHY (NUTS structure)
# ==============================================================================

cat("Assigning geography...\n")

# NUTS1: 5 macro-regions
nuts1_codes <- c("N", "S", "NE", "NO", "C")
nuts1_props <- c(0.25, 0.25, 0.20, 0.15, 0.15)

# Assign NUTS1 to households
households$NUTS1 <- sample(nuts1_codes, nrow(households),
                           replace = TRUE, prob = nuts1_props)

# NUTS2: 6 regions per NUTS1
create_nuts2_code <- function(nuts1) {
  paste0(nuts1, sprintf("%02d", sample(1:6, 1)))
}

households$NUTS2 <- sapply(households$NUTS1, create_nuts2_code)

# NUTS3: 3-4 provinces per NUTS2
create_nuts3_code <- function(nuts2) {
  paste0(nuts2, sprintf("%03d", sample(2:5, 1)))
}

households$NUTS3 <- sapply(households$NUTS2, create_nuts3_code)

# Municipality: 4-8 per NUTS3
create_municipality_code <- function(nuts3) {
  paste0(nuts3, sprintf("%04d", sample(4:12, 1)))
}

households$municipality <- sapply(households$NUTS3, create_municipality_code)

# Merge geography to individuals
individuals <- individuals %>%
  left_join(households %>% select(hh_id, NUTS1, NUTS2, NUTS3, municipality),
            by = "hh_id")

cat("Geography assigned.\n")
cat("NUTS1 distribution:\n")
print(table(households$NUTS1))
cat("\n")

# ==============================================================================
# 8. GENERATE INCOME using model
# ==============================================================================

cat("Generating income...\n")

# Calculate number of employed per household
hh_employment <- individuals %>%
  group_by(hh_id) %>%
  summarise(
    n_employed = sum(employment_status == "Employed", na.rm = TRUE),
    age_head = age[position_in_hh == 1],
    education_head = education_level[position_in_hh == 1],
    .groups = "drop"
  )

# Merge to households
households <- households %>%
  left_join(hh_employment, by = "hh_id")

# Handle missing education
households <- households %>%
  mutate(
    education_head = ifelse(is.na(education_head), "Medium", education_head)
  )

# Generate income using the fitted parameters directly
# Since we fitted the distribution to eq_income, we can use that directly
# and adjust based on household characteristics

# Method 1: If we have the model coefficients, use them
if (!is.null(itsilc_structure$income_model_coef) &&
    all(!is.na(itsilc_structure$income_model_coef))) {

  cat("Using income model coefficients...\n")

  # Create model matrix manually
  households <- households %>%
    mutate(
      # Base prediction (intercept)
      log_eq_income_pred = itsilc_structure$income_model_coef["(Intercept)"],

      # Add education effects (using Medium as reference)
      log_eq_income_pred = log_eq_income_pred +
        ifelse(!is.na(itsilc_structure$income_model_coef["education_headLow"]) & education_head == "Low",
               itsilc_structure$income_model_coef["education_headLow"], 0) +
        ifelse(!is.na(itsilc_structure$income_model_coef["education_headHigh"]) & education_head == "High",
               itsilc_structure$income_model_coef["education_headHigh"], 0),

      # Add other effects
      log_eq_income_pred = log_eq_income_pred +
        ifelse(!is.na(itsilc_structure$income_model_coef["n_employed"]),
               itsilc_structure$income_model_coef["n_employed"] * n_employed, 0) +
        ifelse(!is.na(itsilc_structure$income_model_coef["age_head"]),
               itsilc_structure$income_model_coef["age_head"] * age_head, 0) +
        ifelse(!is.na(itsilc_structure$income_model_coef["I(age_head^2)"]),
               itsilc_structure$income_model_coef["I(age_head^2)"] * age_head^2, 0) +
        ifelse(!is.na(itsilc_structure$income_model_coef["hh_size"]),
               itsilc_structure$income_model_coef["hh_size"] * hh_size, 0),

      # Add residual noise
      log_eq_income = log_eq_income_pred +
        rnorm(n(), 0, itsilc_structure$income_model_sigma),

      # Transform back to original scale
      eq_income = exp(log_eq_income)
    )

} else {

  # Method 2: Use the fitted lognormal distribution with adjustments
  cat("Using direct lognormal generation with adjustments...\n")

  households <- households %>%
    mutate(
      # Base income from fitted distribution
      base_log_income = rnorm(n(),
                              mean = itsilc_structure$income_params["meanlog"],
                              sd = itsilc_structure$income_params["sdlog"]),

      # Adjust for household characteristics
      # Education: High +15%, Low -15%
      edu_adjust = case_when(
        education_head == "High" ~ 0.15,
        education_head == "Low" ~ -0.15,
        TRUE ~ 0
      ),

      # Employment: +10% per employed person
      emp_adjust = n_employed * 0.10,

      # Age: peak at 45-50
      age_adjust = -0.002 * (age_head - 47.5)^2,

      # Household size: small negative effect
      size_adjust = -0.05 * (hh_size - 2),

      # Combine adjustments
      log_eq_income = base_log_income + edu_adjust + emp_adjust + age_adjust + size_adjust,

      # Transform to original scale
      eq_income = exp(log_eq_income)
    )
}

cat("Income generated successfully.\n")

# Recalibrate to match IT-SILC distribution better
# Calculate current log statistics
current_mean <- mean(log(households$eq_income))
current_sd <- sd(log(households$eq_income))

# Target statistics from IT-SILC
target_mean <- itsilc_structure$income_params["meanlog"]
target_sd <- itsilc_structure$income_params["sdlog"]

cat("Recalibrating income distribution...\n")
cat("Current: mean(log)=", round(current_mean, 4), ", sd(log)=", round(current_sd, 4), "\n")
cat("Target:  mean(log)=", round(target_mean, 4), ", sd(log)=", round(target_sd, 4), "\n")

# Rescale in log-space to match target distribution
households <- households %>%
  mutate(
    log_eq_income_rescaled = (log(eq_income) - current_mean) * (target_sd / current_sd) + target_mean,
    eq_income = exp(log_eq_income_rescaled)
  )

cat("After recalibration: mean=", round(mean(households$eq_income)),
    ", median=", round(median(households$eq_income)), "\n")

# Calculate OECD modified equivalence scale
individuals <- individuals %>%
  group_by(hh_id) %>%
  mutate(
    is_first_adult = position_in_hh == 1 & age >= 14,
    is_other_adult = position_in_hh > 1 & age >= 14,
    is_child = age < 14
  ) %>%
  ungroup()

oecd_scales <- individuals %>%
  group_by(hh_id) %>%
  summarise(
    oecd_scale = sum(is_first_adult) * 1.0 +
      sum(is_other_adult) * 0.5 +
      sum(is_child) * 0.3,
    .groups = "drop"
  )

households <- households %>%
  left_join(oecd_scales, by = "hh_id")

# Calculate household total income
households$hh_income <- households$eq_income * households$oecd_scale

# Generate sampling weights
households$weight <- rlnorm(nrow(households),
                            meanlog = itsilc_structure$weight_params["meanlog"],
                            sdlog = itsilc_structure$weight_params["sdlog"])

# Merge income and weights to individuals
individuals <- individuals %>%
  left_join(households %>% select(hh_id, eq_income, hh_income, oecd_scale, weight),
            by = "hh_id")

cat("Income generated.\n")
cat("Equivalent income summary:\n")
print(summary(individuals$eq_income))
cat("\nHousehold income summary:\n")
print(summary(individuals$hh_income))
cat("\n")

# ==============================================================================
# 9. HOUSEHOLD TYPE
# ==============================================================================

individuals <- individuals %>%
  mutate(
    hh_type = case_when(
      hh_size == 1 ~ "Single",
      hh_size == 2 & n_children == 0 ~ "Couple",
      hh_size == 2 & n_children == 1 ~ "Single_parent",
      hh_size >= 3 & n_children == 0 ~ "Other",
      hh_size >= 3 & n_children >= 1 ~ "Family",
      TRUE ~ "Other"
    )
  )

# ==============================================================================
# 10. FINAL DATASET ORGANIZATION
# ==============================================================================

cat("Organizing final dataset...\n")

# Select and reorder columns
synthouse <- individuals %>%
  select(
    # Identifiers
    person_id, hh_id,

    # Geography
    NUTS1, NUTS2, NUTS3, municipality,

    # Individual characteristics
    age, age_class, gender, education_level, employment_status,

    # Household characteristics
    hh_size, hh_type,

    # Income and weights
    eq_income, hh_income, oecd_scale, weight
  ) %>%
  arrange(hh_id, person_id)

# ==============================================================================
# 11. VALIDATION & SUMMARY
# ==============================================================================

cat("\n=== DATASET SUMMARY ===\n")
cat("Total individuals:", nrow(synthouse), "\n")
cat("Total households:", length(unique(synthouse$hh_id)), "\n")
cat("Average household size:", mean(table(synthouse$hh_id)), "\n")
cat("Estimated population size:", sum(unique(synthouse[, c("hh_id", "weight")])$weight), "\n\n")

cat("=== VALIDATION ===\n")
cat("Target avg HH size:", itsilc_structure$mean_hh_size, "\n")
cat("Actual avg HH size:", mean(table(synthouse$hh_id)), "\n\n")

cat("Household size distribution:\n")
print(table(table(synthouse$hh_id)))
cat("\n")

cat("Age class distribution:\n")
print(prop.table(table(synthouse$age_class)))
cat("\n")

cat("Gender distribution:\n")
print(prop.table(table(synthouse$gender)))
cat("\n")

cat("Education distribution (adults):\n")
print(prop.table(table(synthouse$education_level[!is.na(synthouse$education_level)])))
cat("\n")

cat("Employment status:\n")
print(prop.table(table(synthouse$employment_status)))
cat("\n")

cat("NUTS1 distribution:\n")
print(prop.table(table(unique(synthouse[, c("hh_id", "NUTS1")])$NUTS1)))
cat("\n")

# ==============================================================================
# 12. SAVE DATASET
# ==============================================================================

cat("Saving dataset...\n")


# For now: save as RDS
saveRDS(synthouse, "C:/GitHub/inequantiles/synthouse.rds")



# For package development: save as .rda
synthouse <- readRDS("C:/GitHub/inequantiles/synthouse.rds")

# usa compress = "bzip2" per una buona compressione (standard CRAN)
save(synthouse, file = "C:/GitHub/inequantiles/data/synthouse.rda", compress = "bzip2")





