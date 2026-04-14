test_that("qri returns a single numeric scalar", {
  y <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
  result <- qri(y)
  expect_type(result, "double")
  expect_length(result, 1)
})

test_that("qri is 0 for a perfectly equal distribution", {
  # When all values are identical, Q(p/2) = Q(1-p/2) for all p -> QRI = 0
  y <- rep(5, 50)
  expect_equal(qri(y), 0, tolerance = 1e-10)
})

test_that("qri is between 0 and 1 for a positive distribution", {
  set.seed(42)
  y <- rlnorm(200)
  result <- qri(y)
  expect_gte(result, 0)
  expect_lte(result, 1)
})

test_that("qri increases with greater inequality", {
  set.seed(1)
  y_low  <- rlnorm(500, meanlog = 0, sdlog = 0.3)
  y_high <- rlnorm(500, meanlog = 0, sdlog = 1.5)
  expect_lt(qri(y_low), qri(y_high))
})

test_that("qri with equal weights matches unweighted result", {
  set.seed(2)
  y <- rlnorm(100)
  w <- rep(1, 100)
  expect_equal(qri(y, weights = w), qri(y), tolerance = 1e-6)
})

test_that("qri raises a warning for negative values", {
  y <- c(-1, 2, 3, 4, 5)
  expect_warning(qri(y), regexp = "negative")
})

test_that("qri raises an error when length(y) != length(weights)", {
  y <- c(1, 2, 3, 4, 5)
  w <- c(1, 1, 1)
  expect_error(qri(y, weights = w))
})

test_that("qri M parameter: more grid points gives consistent estimate", {
  set.seed(5)
  y <- rlnorm(300)
  r50  <- qri(y, M = 50)
  r200 <- qri(y, M = 200)
  # Both should be close for a smooth distribution
  expect_equal(r50, r200, tolerance = 0.01)
})

test_that("qri na.rm = TRUE handles NAs without error", {
  y <- c(1, 2, NA, 4, 5, 6, 7, 8, 9, 10)
  expect_no_error(qri(y, na.rm = TRUE))
  expect_no_warning(qri(y, na.rm = TRUE))
})
