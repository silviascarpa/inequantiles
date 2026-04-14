test_that("ratio_quantiles returns a single numeric scalar", {
  y <- 1:20
  result <- ratio_quantiles(y)
  expect_type(result, "double")
  expect_length(result, 1)
})

test_that("ratio_quantiles is 1 for a perfectly equal distribution", {
  y <- rep(5, 50)
  expect_equal(as.numeric(ratio_quantiles(y)), 1, tolerance = 1e-10)
})

test_that("ratio_quantiles respects custom prob arguments", {
  set.seed(8)
  y <- rlnorm(200)
  r_90_10 <- ratio_quantiles(y, prob_numerator = 0.90, prob_denominator = 0.10)
  r_75_25 <- ratio_quantiles(y, prob_numerator = 0.75, prob_denominator = 0.25)
  # P90/P10 should be larger than P75/P25 for a right-skewed distribution
  expect_gt(r_90_10, r_75_25)
})

test_that("ratio_quantiles with equal weights matches unweighted result", {
  set.seed(9)
  y <- rlnorm(100)
  w <- rep(1, 100)
  expect_equal(
    ratio_quantiles(y, weights = w),
    ratio_quantiles(y),
    tolerance = 1e-6
  )
})

test_that("ratio_quantiles returns NA with a warning when denominator quantile is 0", {
  y <- c(rep(0, 50), rep(10, 50))
  expect_warning(result <- ratio_quantiles(y, prob_numerator = 0.9, prob_denominator = 0.1))
  expect_true(is.na(result))
})

test_that("ratio_quantiles increases with greater inequality", {
  set.seed(11)
  y_low  <- rlnorm(500, meanlog = 0, sdlog = 0.3)
  y_high <- rlnorm(500, meanlog = 0, sdlog = 1.5)
  expect_lt(ratio_quantiles(y_low), ratio_quantiles(y_high))
})

