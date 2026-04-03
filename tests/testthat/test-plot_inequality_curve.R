test_that("plot_inequality_curve returns a list with p, Rp, qri", {
  result <- plot_inequality_curve(qfunction = qlnorm,
                                  qfun_args = list(meanlog = 9, sdlog = 0.5))
  expect_named(result, c("p", "Rp", "qri"))
})

test_that("plot_inequality_curve p and Rp have length M + 2", {
  result <- plot_inequality_curve(qfunction = qlnorm,
                                  qfun_args = list(meanlog = 9, sdlog = 0.5),
                                  M = 100)
  expect_length(result$p,  102)
  expect_length(result$Rp, 102)
})

test_that("plot_inequality_curve p boundaries are 0 and 1", {
  result <- plot_inequality_curve(qfunction = qlnorm,
                                  qfun_args = list(meanlog = 9, sdlog = 0.5))
  expect_equal(result$p[1],  0)
  expect_equal(result$p[length(result$p)], 1)
})

test_that("plot_inequality_curve Rp boundaries are 0 and 1", {
  result <- plot_inequality_curve(qfunction = qlnorm,
                                  qfun_args = list(meanlog = 9, sdlog = 0.5))
  expect_equal(result$Rp[1],  0)
  expect_equal(result$Rp[length(result$Rp)], 1)
})

test_that("plot_inequality_curve QRI is in (0, 1)", {
  result <- plot_inequality_curve(qfunction = qlnorm,
                                  qfun_args = list(meanlog = 9, sdlog = 0.5))
  expect_gt(result$qri, 0)
  expect_lt(result$qri, 1)
})

test_that("plot_inequality_curve QRI is 0 for a uniform (equal) distribution", {
  # Q(p/2) / Q(1-p/2) = 1 for all p when all values are identical
  result <- plot_inequality_curve(qfunction = qunif,
                                  qfun_args = list(min = 1, max = 1 + 1e-12))
  expect_equal(result$qri, 0, tolerance = 1e-4)
})

test_that("plot_inequality_curve QRI increases with greater inequality", {
  r_low  <- plot_inequality_curve(qfunction = qlnorm,
                                  qfun_args = list(meanlog = 9, sdlog = 0.3))
  r_high <- plot_inequality_curve(qfunction = qlnorm,
                                  qfun_args = list(meanlog = 9, sdlog = 1.2))
  expect_lt(r_low$qri, r_high$qri)
})

test_that("plot_inequality_curve empirical mode returns numeric QRI", {
  set.seed(5)
  y      <- rlnorm(300)
  result <- plot_inequality_curve(y = y)
  expect_true(is.numeric(result$qri))
  expect_length(result$qri, 1)
})

test_that("plot_inequality_curve empirical QRI matches qri()", {
  set.seed(6)
  y      <- rlnorm(300)
  result <- plot_inequality_curve(y = y, M = 100)
  expect_equal(result$qri, round(as.numeric(qri(y, M = 100)), 2))
})

test_that("plot_inequality_curve stops when both y and qfunction provided", {
  expect_error(
    plot_inequality_curve(y = rlnorm(10), qfunction = qlnorm),
    "not both"
  )
})

test_that("plot_inequality_curve stops when neither y nor qfunction provided", {
  expect_error(plot_inequality_curve(), "Provide either")
})

test_that("plot_inequality_curve add = TRUE does not error", {
  plot_inequality_curve(qfunction = qlnorm,
                        qfun_args = list(meanlog = 9, sdlog = 0.5))
  expect_no_error(
    plot_inequality_curve(qfunction = qlnorm,
                          qfun_args = list(meanlog = 9, sdlog = 1.0),
                          col = "tomato", add = TRUE)
  )
})

test_that("plot_inequality_curve returns invisibly", {
  f <- function() plot_inequality_curve(qfunction = qlnorm,
                                        qfun_args = list(meanlog = 9, sdlog = 0.5))
  expect_invisible(f())
})
