test_that("qsr returns a single numeric scalar", {
  y <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
  result <- qsr(y)
  expect_type(result, "double")
  expect_length(result, 1)
})

test_that("qsr is 1 for a perfectly equal distribution", {
  # Top 20% mean == bottom 20% mean when all values identical
  y <- rep(5, 100)
  expect_equal(qsr(y), 1, tolerance = 1e-10)
})

test_that("qsr is >= 1 for a non-negative distribution", {
  set.seed(10)
  y <- rlnorm(200)
  expect_gte(qsr(y), 1)
})

test_that("qsr increases with greater inequality", {
  set.seed(3)
  y_low  <- rlnorm(500, meanlog = 0, sdlog = 0.3)
  y_high <- rlnorm(500, meanlog = 0, sdlog = 1.5)
  expect_lt(qsr(y_low), qsr(y_high))
})

test_that("qsr with equal weights matches unweighted result", {
  set.seed(4)
  y <- rlnorm(100)
  w <- rep(1, 100)
  expect_equal(qsr(y, weights = w), qsr(y), tolerance = 1e-6)
})

test_that("qsr returns NA with a warning when bottom-20% sum is zero", {
  # All bottom values are 0 -> denominator is 0
  y <- c(rep(0, 20), rep(10, 80))
  expect_warning(result <- qsr(y), regexp = "Denominator is zero")
  expect_true(is.na(result))
})

test_that("qsr na.rm = TRUE handles NAs without error", {
  y <- c(1, NA, 3, 4, 5, 6, 7, 8, 9, 10)
  expect_no_error(qsr(y, na.rm = TRUE))
})
