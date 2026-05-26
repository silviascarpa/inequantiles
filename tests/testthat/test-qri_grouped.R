test_that("qri_grouped returns a single numeric scalar", {
  freq  <- c(10, 20, 40, 20, 10)
  lower <- c(0, 10, 20, 30, 40)
  upper <- c(10, 20, 30, 40, 50)
  result <- qri_grouped(freq, lower, upper)
  expect_type(result, "double")
  expect_length(result, 1)
})


test_that("qri_grouped is between 0 and 1 for a typical distribution", {
  freq  <- c(10, 20, 40, 20, 10)
  lower <- c(0, 10, 20, 30, 40)
  upper <- c(10, 20, 30, 40, 50)
  result <- qri_grouped(freq, lower, upper)
  expect_gte(result, 0)
  expect_lte(result, 1)
})

test_that("qri_grouped increases with greater spread", {
  freq_narrow <- c(5, 20, 50, 20, 5)
  freq_wide   <- c(25, 20, 10, 20, 25)
  lower <- c(0, 10, 20, 30, 40)
  upper <- c(10, 20, 30, 40, 50)
  expect_lt(
    qri_grouped(freq_narrow, lower, upper),
    qri_grouped(freq_wide,   lower, upper)
  )
})

test_that("qri_grouped returns NA with a warning for zero total frequency", {
  freq  <- c(0, 0, 0)
  lower <- c(0, 10, 20)
  upper <- c(10, 20, 30)
  expect_warning(
    result <- qri_grouped(freq, lower, upper),
    regexp = "Total frequency is zero"
  )
  expect_true(is.na(result))
})


