test_that("palma_ratio returns a single numeric scalar", {
  y <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
  result <- palma_ratio(y)
  expect_type(result, "double")
  expect_length(result, 1)
})

test_that("palma_ratio is 1 for a perfectly equal distribution", {
  y <- rep(5, 100)
  expect_equal(palma_ratio(y), 1, tolerance = 1e-10)
})

test_that("palma_ratio is >= 1 for a right-skewed distribution", {
  set.seed(20)
  y <- rlnorm(300)
  expect_gte(palma_ratio(y), 1)
})

test_that("palma_ratio increases with greater inequality", {
  set.seed(6)
  y_low  <- rlnorm(500, meanlog = 0, sdlog = 0.3)
  y_high <- rlnorm(500, meanlog = 0, sdlog = 1.5)
  expect_lt(palma_ratio(y_low), palma_ratio(y_high))
})

test_that("palma_ratio with equal weights matches unweighted result", {
  set.seed(7)
  y <- rlnorm(100)
  w <- rep(1, 100)
  expect_equal(palma_ratio(y, weights = w), palma_ratio(y), tolerance = 1e-6)
})

test_that("palma_ratio returns NA with a warning when bottom-40% sum is zero", {
  y <- c(rep(0, 40), rep(10, 60))
  expect_warning(result <- palma_ratio(y))
  expect_true(is.na(result))
})

test_that("palma_ratio na.rm = TRUE handles NAs without error", {
  y <- c(1, NA, 3, 4, 5, 6, 7, 8, 9, 10)
  expect_no_error(palma_ratio(y, na.rm = TRUE))
})
