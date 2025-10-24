# ==============================================================================
# IT-SILC 2024 - Structural Analysis for Synthetic Data Generation
# ==============================================================================
#
# This script analyzes the empirical structure of IT-SILC 2024 to extract
# key patterns and distributions that will be used to generate a realistic
# synthetic dataset for the package.
#
# Input: Itsilc (IT-SILC 2024 dataset)
# Output: List of statistics and models for synthetic data generation
# ==============================================================================

library(dplyr)
library(tidyr)

# Load IT-SILC data (assumed to be in environment as 'Itsilc')
Itsilc <- readRDS("C:/Users/scarp/OneDrive - Unimore/Documenti/Lavori di ricerca/Dati per la ricerca/EU-SILC/ITSILC 24/Dati/Itsilc24.RDS")
Itsilc <- Itsilc[Itsilc$HX050 > 0, ] ## reddito equivalente disponibile
dimHH <- Itsilc %>% group_by(hhid) %>% summarise(nHH = n())

fit_eqInc <- univariateML::mllnorm(Itsilc$HX050)
fit_peso <- univariateML::mllnorm(Itsilc$peso)


cat("IT-SILC 2024 Analysis\n")
cat("=====================\n")
cat("Total individuals:", nrow(Itsilc), "\n")
cat("Total households:", n_distinct(Itsilc$hhid), "\n\n")

# ==============================================================================
# 1. HOUSEHOLD SIZE DISTRIBUTION
# ==============================================================================

cat("1. HOUSEHOLD SIZE DISTRIBUTION\n")
cat("-------------------------------\n")

hh_size_dist <- Itsilc %>%
  group_by(hhid) %>%
  summarise(
    hh_size = n(),
    .groups = "drop"
  ) %>%
  count(hh_size) %>%
  mutate(prop = n / sum(n))

print(hh_size_dist)
cat("Average household size:", mean(rep(hh_size_dist$hh_size, hh_size_dist$n)), "\n\n")

# Save for synthetic generation
hh_size_probs <- setNames(hh_size_dist$prop, hh_size_dist$hh_size)

# ==============================================================================
# 2. AGE DISTRIBUTION
# ==============================================================================

cat("2. AGE DISTRIBUTION\n")
cat("-------------------\n")

# Create age classes
Itsilc <- Itsilc %>%
  mutate(
    age_class = cut(RX010,
                    breaks = c(-1, 14, 17, 24, 34, 49, 64, 100),
                    labels = c("0-14", "15-17", "18-24", "25-34", "35-49", "50-64", "65+"))
  )

age_dist <- table(Itsilc$age_class)
print(prop.table(age_dist))
cat("\n")

# Age distribution by position in household
Itsilc <- Itsilc %>%
  group_by(hhid) %>%
  mutate(position_in_hh = row_number()) %>%
  ungroup()

age_by_position <- Itsilc %>%
  filter(position_in_hh <= 4) %>%  # First 4 members
  group_by(position_in_hh, age_class) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(position_in_hh) %>%
  mutate(prop = n / sum(n))

cat("Age distribution by position in household (first 4 members):\n")
print(age_by_position %>% select(-n) %>% pivot_wider(names_from = age_class, values_from = prop))
cat("\n")

# ==============================================================================
# 3. GENDER DISTRIBUTION
# ==============================================================================

cat("3. GENDER DISTRIBUTION\n")
cat("----------------------\n")

gender_dist <- prop.table(table(Itsilc$RB090))
print(gender_dist)
cat("\n")

# Gender by position and household size
gender_by_position <- Itsilc %>%
  group_by(hhid) %>%
  mutate(hh_size = n()) %>%
  ungroup() %>%
  filter(position_in_hh <= 2, hh_size >= 2) %>%
  group_by(position_in_hh, RB090) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(position_in_hh) %>%
  mutate(prop = n / sum(n))

cat("Gender distribution by position (households with 2+ members):\n")
print(gender_by_position)
cat("\n")

# ==============================================================================
# 4. EDUCATION LEVEL (Macro-categories)
# ==============================================================================

cat("4. EDUCATION LEVEL\n")
cat("------------------\n")

# Create education macro-categories
Itsilc <- Itsilc %>%
  mutate(
    education_macro = case_when(
      PE041 >= 0 & PE041 < 300 ~ "Low",      # No education to lower secondary
      PE041 >= 300 & PE041 < 540 ~ "Medium", # Upper secondary to post-secondary
      PE041 >= 540 ~ "High",                 # Tertiary
      TRUE ~ NA_character_
    )
  )

# Education distribution (adults 18+)
edu_dist <- Itsilc %>%
  filter(RX010 >= 18) %>%
  count(education_macro) %>%
  mutate(prop = n / sum(n))

cat("Education distribution (adults 18+):\n")
print(edu_dist)
cat("\n")

# Education by age class (adults)
edu_by_age <- Itsilc %>%
  filter(RX010 >= 18) %>%
  group_by(age_class, education_macro) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(age_class) %>%
  mutate(prop = n / sum(n)) %>%
  filter(!is.na(education_macro))

cat("Education by age class:\n")
print(edu_by_age %>% select(-n) %>% pivot_wider(names_from = education_macro, values_from = prop))
cat("\n")

# ==============================================================================
# 5. MAIN ACTIVITY STATUS (Occupation)
# ==============================================================================

cat("5. MAIN ACTIVITY STATUS\n")
cat("-----------------------\n")

# Create occupation macro-categories
Itsilc <- Itsilc %>%
  mutate(
    occupation_macro = case_when(
      RB211 == 1 ~ "Employed",
      RB211 == 2 ~ "Unemployed",
      RB211 == 3 ~ "Retired",
      RB211 == 5 ~ "Student",
      RB211 %in% c(4, 6, 7, 8) ~ "Other",
      TRUE ~ NA_character_
    )
  )

occ_dist <- table(Itsilc$occupation_macro)
print(prop.table(occ_dist))
cat("\n")

# Occupation by age class
occ_by_age <- Itsilc %>%
  group_by(age_class, occupation_macro) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(age_class) %>%
  mutate(prop = n / sum(n)) %>%
  filter(!is.na(occupation_macro))

cat("Occupation by age class:\n")
print(occ_by_age %>% select(-n) %>% pivot_wider(names_from = occupation_macro, values_from = prop))
cat("\n")

# Occupation by education (adults 25+)
occ_by_edu <- Itsilc %>%
  filter(RX010 >= 25, !is.na(education_macro)) %>%
  group_by(education_macro, occupation_macro) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(education_macro) %>%
  mutate(prop = n / sum(n)) %>%
  filter(!is.na(occupation_macro))

cat("Occupation by education (adults 25+):\n")
print(occ_by_edu %>% select(-n) %>% pivot_wider(names_from = occupation_macro, values_from = prop))
cat("\n")

# ==============================================================================
# 6. HOUSEHOLD COMPOSITION ARCHETYPES
# ==============================================================================

cat("6. HOUSEHOLD COMPOSITION ARCHETYPES\n")
cat("-----------------------------------\n")

# Define archetypes based on household characteristics
hh_archetypes <- Itsilc %>%
  group_by(hhid) %>%
  summarise(
    hh_size = n(),
    n_children = sum(RX010 < 18),
    n_adults = sum(RX010 >= 18),
    n_elderly = sum(RX010 >= 65),
    age_head = RX010[position_in_hh == 1],
    n_employed = sum(occupation_macro == "Employed", na.rm = TRUE),
    n_students = sum(occupation_macro == "Student", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    archetype = case_when(
      hh_size == 1 & age_head < 35 ~ "Single_young",
      hh_size == 1 & age_head >= 35 & age_head < 65 ~ "Single_adult",
      hh_size == 1 & age_head >= 65 ~ "Single_elderly",
      hh_size == 2 & n_elderly == 2 ~ "Couple_elderly",
      hh_size == 2 & n_children == 0 & age_head < 50 ~ "Couple_young",
      hh_size == 2 & n_children == 0 & age_head >= 50 ~ "Couple_mature",
      hh_size == 2 & n_children == 1 ~ "Parent_child",
      hh_size >= 3 & n_children == 0 ~ "Multi_adult",
      hh_size >= 3 & n_children == 1 ~ "Family_1child",
      hh_size >= 3 & n_children == 2 ~ "Family_2children",
      hh_size >= 3 & n_children >= 3 ~ "Family_3plus_children",
      TRUE ~ "Other"
    )
  )

archetype_dist <- hh_archetypes %>%
  count(archetype) %>%
  mutate(prop = n / sum(n)) %>%
  arrange(desc(n))

cat("Household archetypes (top 15):\n")
print(head(archetype_dist, 15))
cat("\n")

# Save archetype proportions
archetype_probs <- setNames(archetype_dist$prop, archetype_dist$archetype)

# ==============================================================================
# 7. INCOME MODEL
# ==============================================================================

cat("7. INCOME MODEL\n")
cat("---------------\n")

# Prepare data for income model
income_model_data <- Itsilc %>%
  group_by(hhid) %>%
  summarise(
    eq_income = first(HX050),
    peso = first(peso),
    hh_size = n(),
    age_head = RX010[position_in_hh == 1],
    gender_head = RB090[position_in_hh == 1],
    education_head = education_macro[position_in_hh == 1],
    occupation_head = occupation_macro[position_in_hh == 1],
    n_employed = sum(occupation_macro == "Employed", na.rm = TRUE),
    n_children = sum(RX010 < 18),
    .groups = "drop"
  ) %>%
  filter(!is.na(education_head), !is.na(occupation_head))

# Fit income model
income_model <- lm(log(eq_income) ~
                     education_head +
                     n_employed +
                     age_head +
                     I(age_head^2) +
                     hh_size,
                   data = income_model_data)

cat("Income model summary:\n")
print(summary(income_model))
cat("\n")

# Model residual standard deviation
residual_sd <- sigma(income_model)
cat("Residual standard deviation:", residual_sd, "\n\n")

# ==============================================================================
# 8. NUTS1 DISTRIBUTION
# ==============================================================================

cat("8. GEOGRAPHIC DISTRIBUTION (if available)\n")
cat("-----------------------------------------\n")

# Check if NUTS1 or similar geographic variable exists
if ("DB040" %in% colnames(Itsilc)) {
  nuts_dist <- Itsilc %>%
    group_by(hhid) %>%
    summarise(region = first(DB040), .groups = "drop") %>%
    count(region) %>%
    mutate(prop = n / sum(n))

  cat("NUTS1 distribution:\n")
  print(nuts_dist)
} else {
  cat("No geographic variable found. Using default proportions.\n")
}

cat("\n")

# ==============================================================================
# 9. SAVE RESULTS FOR SYNTHETIC DATA GENERATION
# ==============================================================================

cat("9. SAVING RESULTS\n")
cat("-----------------\n")

# Package all results into a list
itsilc_structure <- list(
  # Distributions
  hh_size_probs = hh_size_probs,
  gender_probs = as.vector(gender_dist),
  archetype_probs = archetype_probs,

  # Conditional distributions
  age_by_position = age_by_position,
  gender_by_position = gender_by_position,
  education_by_age = edu_by_age,
  occupation_by_age = occ_by_age,
  occupation_by_education = occ_by_edu,

  # Income model
  income_model_coef = coef(income_model),
  income_model_sigma = residual_sd,

  # Parameters from previous fit
  income_params = c(meanlog = 9.9194, sdlog = 0.6855),
  weight_params = c(meanlog = 6.117978, sdlog = 1.064999),

  # Sample statistics
  mean_hh_size = mean(rep(hh_size_dist$hh_size, hh_size_dist$n)),
  n_households = n_distinct(Itsilc$hhid),
  n_individuals = nrow(Itsilc)
)

# Save to RDS
setwd("C:/GitHub/inequantiles")
saveRDS(itsilc_structure, "itsilc_structure.rds")

cat("Results saved to 'itsilc_structure.rds'\n")
cat("\n=== ANALYSIS COMPLETE ===\n")

# ==============================================================================
# 10. SUMMARY STATISTICS FOR VALIDATION
# ==============================================================================

cat("\n10. SUMMARY FOR VALIDATION\n")
cat("--------------------------\n")

cat("Key statistics to replicate:\n")
cat("- Average household size:", itsilc_structure$mean_hh_size, "\n")
cat("- % Singles:", hh_size_probs["1"] * 100, "%\n")
cat("- % Employed:", mean(Itsilc$occupation_macro == "Employed", na.rm = TRUE) * 100, "%\n")
cat("- % High education (adults):",
    mean(Itsilc$education_macro[Itsilc$RX010 >= 18] == "High", na.rm = TRUE) * 100, "%\n")
cat("- Median eq_income:", median(income_model_data$eq_income), "€\n")
cat("- Mean eq_income:", mean(income_model_data$eq_income), "€\n")
