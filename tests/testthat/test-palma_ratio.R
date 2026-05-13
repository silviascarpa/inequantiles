# Tests for share_ratio() with Palma parameters (prob_numerator=0.90, prob_denominator=0.40)

test_that("share_ratio (Palma) returns a single numeric scalar", {
  y <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
  result <- share_ratio(y, prob_numerator = 0.90, prob_denominator = 0.40)
  expect_type(result, "double")
  expect_length(result, 1)
})

test_that("share_ratio (Palma) is 1 for a perfectly equal distribution", {
  y <- rep(5, 100)
  expect_equal(share_ratio(y, prob_numerator = 0.90, prob_denominator = 0.40),
               1, tolerance = 1e-10)
})

test_that("share_ratio (Palma) is >= 1 for a right-skewed distribution", {
  set.seed(20)
  y <- rlnorm(300)
  expect_gte(share_ratio(y, prob_numerator = 0.90, prob_denominator = 0.40), 1)
})

test_that("share_ratio (Palma) increases with greater inequality", {
  set.seed(6)
  y_low  <- rlnorm(500, meanlog = 0, sdlog = 0.3)
  y_high <- rlnorm(500, meanlog = 0, sdlog = 1.5)
  expect_lt(
    share_ratio(y_low,  prob_numerator = 0.90, prob_denominator = 0.40),
    share_ratio(y_high, prob_numerator = 0.90, prob_denominator = 0.40)
  )
})

test_that("share_ratio (Palma) with equal weights matches unweighted result", {
  set.seed(7)
  y <- rlnorm(100)
  w <- rep(1, 100)
  expect_equal(
    share_ratio(y, weights = w, prob_numerator = 0.90, prob_denominator = 0.40),
    share_ratio(y,              prob_numerator = 0.90, prob_denominator = 0.40),
    tolerance = 1e-6
  )
})

test_that("share_ratio (Palma) returns NA with a warning when bottom-40% sum is zero", {
  y <- c(rep(0, 40), rep(10, 60))
  expect_warning(
    result <- share_ratio(y, prob_numerator = 0.90, prob_denominator = 0.40)
  )
  expect_true(is.na(result))
})

test_that("share_ratio (Palma) na.rm = TRUE handles NAs without error", {
  y <- c(1, NA, 3, 4, 5, 6, 7, 8, 9, 10)
  expect_no_error(share_ratio(y, prob_numerator = 0.90, prob_denominator = 0.40,
                               na.rm = TRUE))
})
