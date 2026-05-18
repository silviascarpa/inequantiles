test_that("share_ratio returns a single numeric scalar", {
  y <- rlnorm(100)
  expect_type(share_ratio(y), "double")
  expect_length(share_ratio(y), 1)
})

test_that("share_ratio is 1 for a perfectly equal distribution", {
  y <- rep(5, 100)
  expect_equal(share_ratio(y), 1, tolerance = 1e-10)
})

test_that("share_ratio increases with greater inequality", {
  set.seed(1)
  y_low  <- rlnorm(500, meanlog = 0, sdlog = 0.3)
  y_high <- rlnorm(500, meanlog = 0, sdlog = 1.5)
  expect_lt(share_ratio(y_low), share_ratio(y_high))
})

test_that("share_ratio with equal weights matches unweighted result", {
  set.seed(2)
  y <- rlnorm(100)
  w <- rep(1, 100)
  expect_equal(share_ratio(y, weights = w), share_ratio(y), tolerance = 1e-6)
})

test_that("share_ratio accepts custom prob_numerator and prob_denominator", {
  y <- rlnorm(100)
  expect_no_error(share_ratio(y, prob_numerator = 0.90, prob_denominator = 0.40))
  expect_no_error(share_ratio(y, prob_numerator = 0.75, prob_denominator = 0.25))
})

test_that("share_ratio returns NA with a warning when denominator is zero", {
  y <- c(rep(0, 20), rep(10, 80))
  expect_warning(result <- share_ratio(y), regexp = "Denominator is zero")
  expect_true(is.na(result))
})

test_that("share_ratio na.rm = TRUE handles NAs without error", {
  y <- c(rlnorm(99), NA)
  expect_no_error(share_ratio(y, na.rm = TRUE))
})
