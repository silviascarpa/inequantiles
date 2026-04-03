test_that("inequantiles returns an object of class 'inequantiles'", {
  y <- rlnorm(100)
  result <- inequantiles(y)
  expect_s3_class(result, "inequantiles")
})

test_that("inequantiles with 'all' returns all four indicators", {
  y <- rlnorm(100)
  result <- inequantiles(y)
  expect_named(result$estimates, c("qri", "qsr", "palma", "p90p10"))
})

test_that("inequantiles subset of indicators returns only requested ones", {
  y <- rlnorm(100)
  result <- inequantiles(y, indicators = c("qri", "palma"))
  expect_named(result$estimates, c("qri", "palma"))
})

test_that("inequantiles single indicator works", {
  y <- rlnorm(100)
  result <- inequantiles(y, indicators = "qsr")
  expect_named(result$estimates, "qsr")
  expect_length(result$estimates, 1)
})

test_that("inequantiles stops on unknown indicator", {
  y <- rlnorm(100)
  expect_error(inequantiles(y, indicators = "gini"), "Unknown indicator")
})

test_that("inequantiles custom prob_num/prob_den changes ratio label", {
  y <- rlnorm(100)
  result <- inequantiles(y, indicators = "ratio_quantiles",
                         prob_num = 0.80, prob_den = 0.20)
  expect_named(result$estimates, "p80p20")
})

test_that("inequantiles estimates are all positive numerics", {
  set.seed(1)
  y <- rlnorm(200)
  result <- inequantiles(y)
  expect_true(all(is.numeric(result$estimates)))
  expect_true(all(result$estimates > 0))
})

test_that("inequantiles se is NULL when se = FALSE", {
  y <- rlnorm(100)
  result <- inequantiles(y, se = FALSE)
  expect_null(result$se)
  expect_null(result$B)
})

test_that("inequantiles stops when se = TRUE but data or strata missing", {
  y <- rlnorm(100)
  expect_error(inequantiles(y, se = TRUE), "'data' must be provided")
  df <- data.frame(y = y, stratum = rep(1:5, 20))
  expect_error(inequantiles(y, se = TRUE, data = df), "'strata' must be provided")
})

test_that("inequantiles with equal weights matches unweighted result", {
  set.seed(2)
  y <- rlnorm(100)
  w <- rep(1, 100)
  r_unweighted <- inequantiles(y)
  r_weighted   <- inequantiles(y, weights = w)
  expect_equal(r_unweighted$estimates, r_weighted$estimates, tolerance = 1e-6)
})

test_that("inequantiles estimates increase with greater inequality", {
  set.seed(3)
  y_low  <- rlnorm(500, meanlog = 0, sdlog = 0.3)
  y_high <- rlnorm(500, meanlog = 0, sdlog = 1.5)
  r_low  <- inequantiles(y_low)
  r_high <- inequantiles(y_high)
  expect_lt(r_low$estimates["qri"],   r_high$estimates["qri"])
  expect_lt(r_low$estimates["qsr"],   r_high$estimates["qsr"])
  expect_lt(r_low$estimates["palma"], r_high$estimates["palma"])
})

test_that("print.inequantiles runs without error", {
  y <- rlnorm(100)
  result <- inequantiles(y)
  expect_output(print(result))
})

test_that("inequantiles se = TRUE returns numeric standard errors", {
  set.seed(4)
  n       <- 200
  strata  <- rep(letters[1:5], each = 40)
  df <- data.frame(y = rlnorm(n), stratum = strata)
  N_h <- c(a = 400, b = 400, c = 400, d = 400, e = 400)
  result <- inequantiles(df$y,
                         se = TRUE, B = 50, seed = 42,
                         data = df, strata = "stratum", N_h = N_h)
  expect_length(result$se, 4)
  expect_true(all(is.numeric(result$se)))
  expect_true(all(result$se >= 0))
  expect_equal(result$B, 50)
})
