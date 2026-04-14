test_that("quantile_grouped returns a named numeric vector", {
  freq   <- c(10, 20, 30, 25, 15)
  lower  <- c(0, 10, 20, 30, 40)
  upper  <- c(10, 20, 30, 40, 50)
  result <- quantile_grouped(freq, lower, upper, probs = c(0.25, 0.5, 0.75))
  expect_type(result, "double")
  expect_length(result, 3)
  expect_named(result)
})

test_that("quantile_grouped median lies within the correct class interval", {
  freq  <- c(0, 0, 100, 0, 0)
  lower <- c(0, 10, 20, 30, 40)
  upper <- c(10, 20, 30, 40, 50)
  result <- quantile_grouped(freq, lower, upper, probs = 0.5)
  expect_gte(result, 20)
  expect_lte(result, 30)
})

test_that("quantile_grouped handles open lower bound (-Inf) without error", {
  freq  <- c(10, 20, 30, 25, 15)
  lower <- c(-Inf, 10, 20, 30, 40)
  upper <- c(10, 20, 30, 40, 50)
  expect_no_error(quantile_grouped(freq, lower, upper, probs = 0.5))
})

test_that("quantile_grouped handles open upper bound (Inf) without error", {
  freq  <- c(10, 20, 30, 25, 15)
  lower <- c(0, 10, 20, 30, 40)
  upper <- c(10, 20, 30, 40, Inf)
  expect_no_error(quantile_grouped(freq, lower, upper, probs = 0.5))
})

test_that("quantile_grouped returns NA for zero total frequency", {
  freq  <- c(0, 0, 0)
  lower <- c(0, 10, 20)
  upper <- c(10, 20, 30)
  result <- quantile_grouped(freq, lower, upper, probs = 0.5)
  expect_true(is.na(result))
})

test_that("quantile_grouped errors when vectors have different lengths", {
  expect_error(
    quantile_grouped(c(10, 20), c(0, 10, 20), c(10, 20, 30), probs = 0.5)
  )
})

test_that("quantile_grouped result is monotone in probs", {
  freq  <- c(5, 15, 40, 30, 10)
  lower <- c(0, 10, 20, 30, 40)
  upper <- c(10, 20, 30, 40, 50)
  result <- quantile_grouped(freq, lower, upper, probs = c(0.1, 0.25, 0.5, 0.75, 0.9))
  expect_true(all(diff(result) >= 0))
})
