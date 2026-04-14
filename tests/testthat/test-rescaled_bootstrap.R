# Helper: build a minimal stratified dataset for simple (SRS) design
make_srs_data <- function(seed = 1) {
  set.seed(seed)
  n <- 60
  data.frame(
    income  = rlnorm(n, meanlog = 3, sdlog = 0.5),
    stratum = rep(c("A", "B", "C"), each = 20)
  )
}

# Helper: build a minimal dataset for complex (PSU-weighted) design
make_psu_data <- function(seed = 2) {
  set.seed(seed)
  n <- 60
  data.frame(
    income  = rlnorm(n, meanlog = 3, sdlog = 0.5),
    stratum = rep(c("A", "B", "C"), each = 20),
    psu     = rep(1:6, each = 10),
    weight  = runif(n, 0.5, 2)
  )
}

# --------------------------------------------------------------------------
# Simple design (SRS)
# --------------------------------------------------------------------------

test_that("rescaled_bootstrap (simple) returns correct class and components", {
  dat <- make_srs_data()
  result <- rescaled_bootstrap(
    data      = dat,
    y         = "income",
    strata    = "stratum",
    N_h       = c(A = 200, B = 200, C = 200),
    estimator = function(y) mean(y),
    B         = 50,
    seed      = 1
  )
  expect_s3_class(result, "rescaled_bootstrap")
  expect_true(all(c("variance", "boot_estimates", "B", "design") %in% names(result)))
})

test_that("rescaled_bootstrap (simple, by_strata = FALSE) variance is a positive scalar", {
  dat <- make_srs_data()
  result <- rescaled_bootstrap(
    data      = dat,
    y         = "income",
    strata    = "stratum",
    N_h       = c(A = 200, B = 200, C = 200),
    estimator = function(y) mean(y),
    B         = 50,
    seed      = 1,
    by_strata = FALSE
  )
  expect_type(result$variance, "double")
  expect_length(result$variance, 1)
  expect_gt(result$variance, 0)
})

test_that("rescaled_bootstrap (simple, by_strata = TRUE) variance is a positive vector", {
  dat <- make_srs_data()
  result <- rescaled_bootstrap(
    data      = dat,
    y         = "income",
    strata    = "stratum",
    N_h       = c(A = 200, B = 200, C = 200),
    estimator = function(y) mean(y),
    B         = 50,
    seed      = 1,
    by_strata = TRUE
  )
  expect_type(result$variance, "double")
  expect_true(all(result$variance > 0))
})

test_that("rescaled_bootstrap (simple) B replicate estimates have correct length", {
  dat <- make_srs_data()
  B <- 50
  result <- rescaled_bootstrap(
    data      = dat,
    y         = "income",
    strata    = "stratum",
    N_h       = c(A = 200, B = 200, C = 200),
    estimator = function(y) mean(y),
    B         = B,
    seed      = 1,
    by_strata = FALSE
  )
  expect_length(result$boot_estimates, B)
})

test_that("rescaled_bootstrap (simple) errors when N_h is missing", {
  dat <- make_srs_data()
  expect_error(
    rescaled_bootstrap(
      data      = dat,
      y         = "income",
      strata    = "stratum",
      estimator = function(y) mean(y),
      B         = 20
    )
  )
})

# --------------------------------------------------------------------------
# Complex design (PSU + weights)
# --------------------------------------------------------------------------

test_that("rescaled_bootstrap (complex) returns correct class and components", {
  dat <- make_psu_data()
  result <- rescaled_bootstrap(
    data      = dat,
    y         = "income",
    strata    = "stratum",
    psu       = "psu",
    weights   = "weight",
    estimator = function(y, w) weighted.mean(y, w),
    B         = 50,
    seed      = 2
  )
  expect_s3_class(result, "rescaled_bootstrap")
  expect_true(all(c("variance", "boot_estimates", "B", "design") %in% names(result)))
})

test_that("rescaled_bootstrap (complex, by_strata = FALSE) variance is a positive scalar", {
  dat <- make_psu_data()
  result <- rescaled_bootstrap(
    data      = dat,
    y         = "income",
    strata    = "stratum",
    psu       = "psu",
    weights   = "weight",
    estimator = function(y, w) weighted.mean(y, w),
    B         = 50,
    seed      = 2,
    by_strata = FALSE
  )
  expect_length(result$variance, 1)
  expect_gt(result$variance, 0)
})

test_that("rescaled_bootstrap (complex, by_strata = TRUE) variance is a positive vector", {
  dat <- make_psu_data()
  result <- rescaled_bootstrap(
    data      = dat,
    y         = "income",
    strata    = "stratum",
    psu       = "psu",
    weights   = "weight",
    estimator = function(y, w) weighted.mean(y, w),
    B         = 50,
    seed      = 2,
    by_strata = TRUE
  )
  expect_type(result$variance, "double")
  expect_true(all(result$variance > 0))
})

test_that("rescaled_bootstrap errors when weights given without psu", {
  dat <- make_psu_data()
  expect_error(
    rescaled_bootstrap(
      data      = dat,
      y         = "income",
      strata    = "stratum",
      weights   = "weight",   # no psu
      estimator = function(y, w) weighted.mean(y, w),
      B         = 20
    )
  )
})

test_that("rescaled_bootstrap errors when psu given without weights", {
  dat <- make_psu_data()
  expect_error(
    rescaled_bootstrap(
      data      = dat,
      y         = "income",
      strata    = "stratum",
      psu       = "psu",      # no weights
      estimator = function(y, w) weighted.mean(y, w),
      B         = 20
    )
  )
})

test_that("rescaled_bootstrap seed argument makes results reproducible", {
  dat <- make_srs_data()
  args <- list(
    data      = dat,
    y         = "income",
    strata    = "stratum",
    N_h       = c(A = 200, B = 200, C = 200),
    estimator = function(y) mean(y),
    B         = 30,
    seed      = 99,
    by_strata = FALSE
  )
  r1 <- do.call(rescaled_bootstrap, args)
  r2 <- do.call(rescaled_bootstrap, args)
  expect_equal(r1$variance, r2$variance)
})

