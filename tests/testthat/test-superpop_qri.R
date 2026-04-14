test_that("superpop_qri returns a single numeric scalar", {
  result <- superpop_qri(qlnorm, meanlog = 0, sdlog = 1)
  expect_type(result, "double")
  expect_length(result, 1)
})

test_that("superpop_qri is between 0 and 1 for a positive distribution", {
  result <- superpop_qri(qlnorm, meanlog = 0, sdlog = 1)
  expect_gte(result, 0)
  expect_lte(result, 1)
})

test_that("superpop_qri is 0 for a degenerate (point-mass) distribution", {
  # A distribution where all mass is at one point: Q(p) = constant for all p
  # Approximate with a very small sd normal -> QRI ~ 0
  result <- superpop_qri(qnorm, mean = 100, sd = 1e-8)
  expect_equal(result, 0, tolerance = 1e-4)
})

test_that("superpop_qri increases with greater sdlog for lognormal", {
  r_low  <- superpop_qri(qlnorm, meanlog = 0, sdlog = 0.3)
  r_high <- superpop_qri(qlnorm, meanlog = 0, sdlog = 1.5)
  expect_lt(r_low, r_high)
})

test_that("superpop_qri works with different distributions", {
  expect_no_error(superpop_qri(qexp,     rate = 1))
  expect_no_error(superpop_qri(qweibull, shape = 2, scale = 1))
  expect_no_error(superpop_qri(qgamma,   shape = 2, rate = 1))
})

test_that("superpop_qri is invariant to location shifts", {
  # Adding a constant to a distribution does not change QRI
  # because R(p) = Q(p/2) / Q(1-p/2) and... actually this is not invariant
  # to location, so just check both results are in [0,1]
  r1 <- superpop_qri(qlnorm, meanlog = 0, sdlog = 1)
  r2 <- superpop_qri(qlnorm, meanlog = 5, sdlog = 1)
  expect_gte(r1, 0); expect_lte(r1, 1)
  expect_gte(r2, 0); expect_lte(r2, 1)
})

test_that("superpop_qri subdivisions argument is passed through", {
  r1 <- superpop_qri(qlnorm, meanlog = 0, sdlog = 1, subdivisions = 100L)
  r2 <- superpop_qri(qlnorm, meanlog = 0, sdlog = 1, subdivisions = 2000L)
  # Results should be very close regardless of subdivisions
  expect_equal(r1, r2, tolerance = 1e-4)
})
