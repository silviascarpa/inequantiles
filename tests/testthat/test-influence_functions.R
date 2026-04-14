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
# if_qsr
# --------------------------------------------------------------------------

test_that("if_qsr returns a numeric vector of the same length as y", {
  result <- if_qsr(y_small)
  expect_type(result, "double")
  expect_length(result, length(y_small))
})

test_that("if_qsr (weighted) returns finite values of correct length", {
  result <- if_qsr(y_small, weights = w_small)
  expect_length(result, length(y_small))
  expect_true(all(is.finite(result)))
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
