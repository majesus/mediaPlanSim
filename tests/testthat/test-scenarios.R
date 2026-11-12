r1 <- c(A = 0.20, B = 0.15)
r2 <- c(A = 0.312, B = 0.23925)

test_that("constructors return classed scenarios with names and insertions", {
  s <- scenario_independent(r1, r2, c(5, 4))
  expect_s3_class(s, "mp_scenario")
  expect_equal(s$names, c("A", "B"))
  expect_equal(s$insertions, c(5L, 4L))
  s2 <- scenario_gaussian(unname(r1), unname(r2), c(5, 4), rho = 0.3)
  expect_equal(s2$names, c("V1", "V2"))
  expect_equal(unname(s2$rho[1, 2]), 0.3)
  s3 <- scenario_latent_class(c(2, 2), matrix(c(.4, .1, .1, .4), 2), c(3, 3))
  expect_equal(sum(s3$class_prob), 1)
  expect_equal(s3$names, c("V1", "V2"))
})

test_that("invalid marginals are rejected with informative errors", {
  expect_error(scenario_independent(c(A = 0.2), c(A = 0.5), 5), "r2 must lie")
  expect_error(scenario_independent(c(A = 0.2), c(A = 0.2), 5), "r2 must lie")
  expect_error(scenario_independent(c(A = 20), c(A = 31.2), 5), "proportions")
  expect_error(scenario_independent(r1, r2[1], c(5, 5)), "equal length")
  expect_error(scenario_independent(r1, r2, c(5, 2.5)), "positive integers")
  expect_error(scenario_independent(r1, r2, c(0, 5)), "positive integers")
  expect_error(scenario_independent(r1, r2, c(300, 200)), "limited to 400")
  expect_error(scenario_independent(c(A = 0.2, A = 0.1), c(0.3, 0.19), c(5, 5)),
               "unique")
  expect_error(scenario_independent("a", r2, c(5, 5)), "numeric")
})

test_that("rho is validated", {
  expect_error(scenario_gaussian(r1, r2, c(5, 5), rho = 2), "valid")
  expect_error(scenario_gaussian(r1, r2, c(5, 5), rho = "a"), "single number")
  bad <- matrix(c(1, 0.9, 0.9, 1), 2)
  bad2 <- matrix(c(1, 0.9, -0.9, 1), 2)
  expect_s3_class(scenario_gaussian(r1, r2, c(5, 5), rho = bad), "mp_scenario")
  expect_error(scenario_gaussian(r1, r2, c(5, 5), rho = bad2), "symmetric")
  expect_error(scenario_gaussian(r1, r2, c(5, 5), rho = diag(3)), "m x m")
  r3 <- matrix(c(1, .9, -.9, .9, 1, .9, -.9, .9, 1), 3)
  expect_error(scenario_gaussian(c(.2, .15, .1), c(.312, .23925, .17), c(5, 5, 5),
                                 rho = r3), "positive semi-definite")
  expect_error(scenario_gaussian(c(.2, .15, .1), c(.312, .23925, .17), c(5, 5, 5),
                                 rho = -0.9), "positive semi-definite")
})

test_that("latent-class inputs are validated", {
  p <- matrix(c(.4, .1, .1, .4), 2)
  expect_error(scenario_latent_class(c(.5, .5, 0), p, c(3, 3)), "one row")
  expect_error(scenario_latent_class(c(.5, -.5), p, c(3, 3)), "non-negative")
  expect_error(scenario_latent_class(c(.5, .5), p * 3, c(3, 3)), "\\[0, 1\\]")
  expect_error(scenario_latent_class(c(.5, .5), p, c(3, 3, 3)), "positive integers")
  expect_error(scenario_latent_class(c(.5, NA), p, c(3, 3)), "NAs")
})

test_that("scenarios print", {
  expect_output(print(scenario_gaussian(r1, r2, c(5, 5), rho = 0.2)), "gaussian")
  expect_output(print(scenario_latent_class(c(.5, .5), matrix(c(.4, .1, .1, .4), 2),
                                            c(3, 3))), "Segments")
})
