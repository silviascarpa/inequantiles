# Shared test data
set.seed(42)
y_small <- rlnorm(80, meanlog = 3, sdlog = 0.5)
w_small <- runif(80, 0.5, 2)

# --------------------------------------------------------------------------
# if_quantile
# --------------------------------------------------------------------------

test_that("if_quantile returns a numeric vector of the same length as y", {
  result <- if_quantile(y_small, probs = 0.5)
  expect_type(result, "double")
  expect_length(result, length(y_small))
})

test_that("if_quantile (weighted) returns finite values of correct length", {
  result <- if_quantile(y_small, weights = w_small, probs = 0.5)
  expect_length(result, length(y_small))
  expect_true(all(is.finite(result)))
})

test_that("if_quantile accepts different prob values", {
  expect_no_error(if_quantile(y_small, probs = 0.25))
  expect_no_error(if_quantile(y_small, probs = 0.75))
  expect_no_error(if_quantile(y_small, probs = 0.90))
})

test_that("if_quantile na.rm = TRUE handles NAs", {
  y_na <- c(y_small, NA)
  expect_no_error(if_quantile(y_na, probs = 0.5, na.rm = TRUE))
})

# --------------------------------------------------------------------------
# if_qri
# --------------------------------------------------------------------------

test_that("if_qri returns a numeric vector of the same length as y", {
  result <- if_qri(y_small)
  expect_type(result, "double")
  expect_length(result, length(y_small))
})

test_that("if_qri (weighted) returns finite values of correct length", {
  result <- if_qri(y_small, weights = w_small)
  expect_length(result, length(y_small))
  expect_true(all(is.finite(result)))
})

# --------------------------------------------------------------------------
# if_share_ratio
# --------------------------------------------------------------------------

test_that("if_share_ratio returns a numeric vector of the same length as y", {
  result <- if_share_ratio(y_small)
  expect_type(result, "double")
  expect_length(result, length(y_small))
})

test_that("if_share_ratio (weighted) returns finite values of correct length", {
  result <- if_share_ratio(y_small, weights = w_small)
  expect_length(result, length(y_small))
  expect_true(all(is.finite(result)))
})

test_that("if_share_ratio accepts custom prob_numerator and prob_denominator", {
  expect_no_error(if_share_ratio(y_small, prob_numerator = 0.90,
                                 prob_denominator = 0.40))
  expect_no_error(if_share_ratio(y_small, prob_numerator = 0.75,
                                 prob_denominator = 0.25))
})

test_that("if_share_ratio na.rm = TRUE handles NAs without error", {
  y_na <- c(y_small, NA)
  expect_no_error(if_share_ratio(y_na, na.rm = TRUE))
  expect_length(if_share_ratio(y_na, na.rm = TRUE), length(y_small))
})

# --------------------------------------------------------------------------
# if_gini
# --------------------------------------------------------------------------

test_that("if_gini returns a numeric vector of the same length as y", {
  result <- if_gini(y_small)
  expect_type(result, "double")
  expect_length(result, length(y_small))
})

test_that("if_gini (weighted) returns finite values of correct length", {
  result <- if_gini(y_small, weights = w_small)
  expect_length(result, length(y_small))
  expect_true(all(is.finite(result)))
})

test_that("if_gini errors when weights contain negative values", {
  w_neg <- w_small
  w_neg[1] <- -1
  expect_error(if_gini(y_small, weights = w_neg))
})

test_that("if_gini errors when length(weights) != length(y)", {
  # suppressWarnings() silences the recycling warning that R emits before
  # if_gini reaches its own length check and stops with an error.
  expect_error(suppressWarnings(if_gini(y_small, weights = w_small[-1])))
})

test_that("if_gini na.rm = TRUE handles NAs without error", {
  y_na <- c(y_small, NA)
  expect_no_error(if_gini(y_na, na.rm = TRUE))
})

# --------------------------------------------------------------------------
# if_ratio_quantiles
# --------------------------------------------------------------------------

test_that("if_ratio_quantiles returns a numeric vector of the same length as y", {
  result <- if_ratio_quantiles(y_small)
  expect_type(result, "double")
  expect_length(result, length(y_small))
})

test_that("if_ratio_quantiles (weighted) returns finite values of correct length", {
  result <- if_ratio_quantiles(y_small, weights = w_small)
  expect_length(result, length(y_small))
  expect_true(all(is.finite(result)))
})

test_that("if_ratio_quantiles accepts custom prob_numerator and prob_denominator", {
  expect_no_error(if_ratio_quantiles(y_small, prob_numerator = 0.80,
                                     prob_denominator = 0.20))
  expect_no_error(if_ratio_quantiles(y_small, prob_numerator = 0.75,
                                     prob_denominator = 0.25))
})

test_that("if_ratio_quantiles sums to approximately zero (mean-zero property)", {
  result <- if_ratio_quantiles(y_small)
  expect_lt(abs(mean(result)), 0.01)
})

test_that("if_ratio_quantiles (weighted) sums to approximately zero", {
  result <- if_ratio_quantiles(y_small, weights = w_small)
  expect_lt(abs(weighted.mean(result, w_small)), 0.01)
})

test_that("if_ratio_quantiles na.rm = TRUE handles NAs without error", {
  y_na <- c(y_small, NA)
  expect_no_error(if_ratio_quantiles(y_na, na.rm = TRUE))
  expect_length(if_ratio_quantiles(y_na, na.rm = TRUE), length(y_small))
})

