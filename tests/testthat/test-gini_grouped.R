test_that("gini_grouped returns a single numeric scalar", {
  Y    <- c(100, 200, 300, 400, 500)
  freq <- c(20, 20, 20, 20, 20)
  result <- gini_grouped(Y, freq)
  expect_type(result, "double")
  expect_length(result, 1)
})

test_that("gini_grouped is between 0 and 1", {
  Y    <- c(100, 200, 400, 800)
  freq <- c(40, 30, 20, 10)
  result <- gini_grouped(Y, freq)
  expect_gte(result, 0)
  expect_lte(result, 1)
})

test_that("gini_grouped is 0 for a perfectly equal distribution", {
  # Equal Y and equal freq -> perfect equality
  Y    <- c(100, 100, 100, 100)
  freq <- c(25, 25, 25, 25)
  result <- gini_grouped(Y, freq)
  expect_equal(result, 0, tolerance = 1e-10)
})

test_that("gini_grouped increases with greater inequality", {
  freq <- c(25, 25, 25, 25)
  Y_low  <- c(80, 90, 110, 120)
  Y_high <- c(10, 20, 100, 870)
  expect_lt(gini_grouped(Y_low, freq), gini_grouped(Y_high, freq))
})

test_that("gini_grouped errors when Y and freq have different lengths", {
  expect_error(gini_grouped(c(100, 200, 300), c(10, 20)))
})

test_that("gini_grouped errors when any freq is negative", {
  expect_error(gini_grouped(c(100, 200), c(10, -5)))
})

test_that("gini_grouped errors when total frequency is 0", {
  expect_error(gini_grouped(c(100, 200), c(0, 0)))
})

test_that("gini_grouped warns (does not stop) for negative Y values", {
  Y    <- c(-50, 100, 200, 300)
  freq <- c(10, 30, 30, 30)
  expect_warning(gini_grouped(Y, freq))
})

test_that("gini_grouped returns 0 with a warning when sum(Y) is 0", {
  Y    <- c(0, 0, 0)
  freq <- c(10, 20, 30)
  expect_warning(result <- gini_grouped(Y, freq))
  expect_equal(result, 0)
})
