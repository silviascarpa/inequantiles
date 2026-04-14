test_that("csquantile returns a named numeric vector", {
  y <- c(1, 2, 3, 4, 5)
  result <- csquantile(y, probs = c(0.25, 0.5, 0.75))
  expect_type(result, "double")
  expect_length(result, 3)
  expect_named(result)
})

test_that("csquantile without weights matches base R quantile for types 4-9", {
  set.seed(42)
  y <- rnorm(100)
  probs <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  for (tp in 4:9) {
    expect_equal(
      unname(csquantile(y, probs = probs, type = tp)),
      unname(quantile(y, probs = probs, type = tp)),
      tolerance = 1e-10,
      label = paste("type =", tp)
    )
  }
})

test_that("csquantile with equal weights matches unweighted result", {
  set.seed(1)
  y <- rnorm(50)
  w <- rep(1, 50)
  probs <- c(0.25, 0.5, 0.75)
  expect_equal(
    unname(csquantile(y, weights = w, probs = probs, type = 6)),
    unname(csquantile(y, weights = NULL, probs = probs, type = 6)),
    tolerance = 1e-6
  )
})

test_that("csquantile boundary: p = 0 returns minimum, p ~ 1 returns maximum", {
  y <- c(3, 1, 4, 1, 5, 9, 2, 6)
  expect_equal(unname(csquantile(y, probs = 0)), min(y))
  expect_equal(unname(csquantile(y, probs = 1)), max(y))
})

test_that("csquantile type HD returns a value in [min(y), max(y)]", {
  set.seed(7)
  y <- rexp(80)
  result <- csquantile(y, probs = c(0.1, 0.5, 0.9), type = "HD")
  expect_true(all(result >= min(y)))
  expect_true(all(result <= max(y)))
  expect_length(result, 3)
})

test_that("csquantile type HD with weights returns a value in [min(y), max(y)]", {
  set.seed(7)
  y <- rexp(50)
  w <- runif(50, 0.5, 2)
  result <- csquantile(y, weights = w, probs = c(0.25, 0.75), type = "HD")
  expect_true(all(result >= min(y)))
  expect_true(all(result <= max(y)))
})

test_that("csquantile na.rm = TRUE removes NAs without error", {
  y <- c(1, 2, NA, 4, 5)
  expect_no_error(csquantile(y, probs = 0.5, na.rm = TRUE))
})

test_that("csquantile na.rm = FALSE errors when NAs are present (consistent with base R)", {
  y <- c(1, 2, NA, 4, 5)
  expect_error(csquantile(y, probs = 0.5, na.rm = FALSE))
})

test_that("csquantile result is monotone in probs", {
  set.seed(3)
  y <- rnorm(200)
  probs <- seq(0.1, 0.9, by = 0.1)
  result <- csquantile(y, probs = probs, type = 7)
  expect_true(all(diff(result) >= 0))
})

test_that("csquantile weighted quantile: heavier weights pull estimate toward high values", {
  y <- c(1, 10)
  # weight all on 10 -> median should be 10
  w <- c(0.001, 999)
  expect_gt(csquantile(y, weights = w, probs = 0.5, type = 4), 5)
})

