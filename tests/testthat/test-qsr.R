# Tests for share_ratio() with QSR parameters (prob_numerator=0.80, prob_denominator=0.20)

test_that("share_ratio (QSR) returns a single numeric scalar", {
  y <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
  result <- share_ratio(y, prob_numerator = 0.80, prob_denominator = 0.20)
  expect_type(result, "double")
  expect_length(result, 1)
})

test_that("share_ratio (QSR) is 1 for a perfectly equal distribution", {
  y <- rep(5, 100)
  expect_equal(share_ratio(y, prob_numerator = 0.80, prob_denominator = 0.20),
               1, tolerance = 1e-10)
})

test_that("share_ratio (QSR) is >= 1 for a non-negative distribution", {
  set.seed(10)
  y <- rlnorm(200)
  expect_gte(share_ratio(y, prob_numerator = 0.80, prob_denominator = 0.20), 1)
})

test_that("share_ratio (QSR) increases with greater inequality", {
  set.seed(3)
  y_low  <- rlnorm(500, meanlog = 0, sdlog = 0.3)
  y_high <- rlnorm(500, meanlog = 0, sdlog = 1.5)
  expect_lt(
    share_ratio(y_low,  prob_numerator = 0.80, prob_denominator = 0.20),
    share_ratio(y_high, prob_numerator = 0.80, prob_denominator = 0.20)
  )
})

test_that("share_ratio (QSR) with equal weights matches unweighted result", {
  set.seed(4)
  y <- rlnorm(100)
  w <- rep(1, 100)
  expect_equal(
    share_ratio(y, weights = w, prob_numerator = 0.80, prob_denominator = 0.20),
    share_ratio(y,              prob_numerator = 0.80, prob_denominator = 0.20),
    tolerance = 1e-6
  )
})

test_that("share_ratio (QSR) returns NA with a warning when bottom-20% sum is zero", {
  y <- c(rep(0, 20), rep(10, 80))
  expect_warning(
    result <- share_ratio(y, prob_numerator = 0.80, prob_denominator = 0.20),
    regexp = "Denominator is zero"
  )
  expect_true(is.na(result))
})

test_that("share_ratio (QSR) na.rm = TRUE handles NAs without error", {
  y <- c(1, NA, 3, 4, 5, 6, 7, 8, 9, 10)
  expect_no_error(share_ratio(y, prob_numerator = 0.80, prob_denominator = 0.20,
                               na.rm = TRUE))
})
