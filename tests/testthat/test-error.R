test_that("distribution_error is zero for equal distributions", {
  p <- c(0.3, 0.4, 0.2, 0.1)
  e <- distribution_error(p, p)
  expect_equal(unname(e), rep(0, 5))
  expect_equal(names(e), c("tv", "ks", "reach_1plus", "reach_3plus", "mean_contacts"))
})

test_that("distribution_error computes the documented quantities", {
  a <- c(0.5, 0.3, 0.2)
  b <- c(0.4, 0.4, 0.1, 0.1)
  e <- distribution_error(a, b)
  expect_equal(e[["tv"]], 0.5 * (0.1 + 0.1 + 0.1 + 0.1))
  expect_equal(e[["ks"]], 0.1)
  expect_equal(e[["reach_1plus"]], 0.5 - 0.6)
  expect_equal(e[["reach_3plus"]], (1 - 1) - (1 - 0.9))
  expect_equal(e[["mean_contacts"]], (0.3 + 0.4) - (0.4 + 0.2 + 0.3))
})

test_that("inputs of several kinds are accepted and invalid ones rejected", {
  d <- data.frame(k = 0:2, pmf = c(0.5, 0.3, 0.2))
  expect_equal(distribution_error(d, list(distribution = d))[["tv"]], 0)
  expect_error(distribution_error(c(0.5, 0.6), c(0.5, 0.5)), "sum to one")
  expect_error(distribution_error("a", c(0.5, 0.5)), "numeric")
  expect_error(distribution_error(c(1.2, -0.2), c(0.5, 0.5)), "numeric")
})
