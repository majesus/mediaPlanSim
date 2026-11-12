r1 <- c(A = 0.20, B = 0.15)
r2 <- c(A = 0.312, B = 0.23925)
scns <- list(
  independent = scenario_independent(r1, r2, c(5, 4)),
  gaussian = scenario_gaussian(r1, r2, c(5, 4), rho = 0.5),
  latent = scenario_latent_class(c(.35, .35, .3),
                                 cbind(A = c(.45, .15, .15), B = c(.15, .45, .15)),
                                 c(5, 4))
)

test_that("simulation has the right shape, support and names", {
  for (s in scns) {
    sim <- simulate_exposure(s, 500, seed = 1)
    expect_s3_class(sim, "mp_sim")
    expect_equal(dim(sim$contacts), c(500L, 2L))
    expect_equal(colnames(sim$contacts), c("A", "B"))
    expect_true(all(sim$contacts[, 1] >= 0 & sim$contacts[, 1] <= s$insertions[1]))
    expect_true(all(sim$contacts[, 2] >= 0 & sim$contacts[, 2] <= s$insertions[2]))
  }
})

test_that("seeds give reproducible draws and restore the RNG state", {
  s <- scns$gaussian
  a <- simulate_exposure(s, 300, seed = 42)
  b <- simulate_exposure(s, 300, seed = 42)
  expect_identical(a$contacts, b$contacts)
  expect_false(identical(a$contacts, simulate_exposure(s, 300, seed = 43)$contacts))
  set.seed(7)
  x <- runif(1)
  set.seed(7)
  invisible(simulate_exposure(s, 300, seed = 99))
  expect_equal(runif(1), x)
})

test_that("sample moments converge to the population moments", {
  for (s in scns) {
    sim <- simulate_exposure(s, 300000, seed = 5)
    sm <- sample_moments(sim)
    pm <- population_moments(s)
    expect_lt(max(abs(sm$r1 - pm$r1)), 3e-3)
    expect_lt(max(abs(sm$r2 - pm$r2)), 4e-3)
    expect_lt(abs(sm$dup[1, 2] - pm$dup[1, 2]), 3e-3)
  }
})

test_that("the empirical distribution converges to the truth", {
  for (s in scns) {
    sim <- simulate_exposure(s, 200000, seed = 8)
    err <- distribution_error(empirical_distribution(sim), exposure_truth(s))
    expect_lt(err[["tv"]], 0.01)
    expect_lt(abs(err[["reach_1plus"]]), 0.006)
  }
})

test_that("r2 is NA for a vehicle with a single insertion", {
  s <- scenario_latent_class(c(.5, .5), cbind(A = c(.4, .1), B = c(.1, .4)), c(1, 3))
  sm <- sample_moments(simulate_exposure(s, 1000, seed = 1))
  expect_true(is.na(sm$r2["A"]))
  expect_false(is.na(sm$r2["B"]))
})

test_that("export and printing work", {
  sim <- simulate_exposure(scns$independent, 50, seed = 1)
  df <- as.data.frame(sim)
  expect_equal(names(df), c("id", "A", "B", "total"))
  expect_equal(df$total, rowSums(sim$contacts))
  expect_output(print(sim), "50 people")
  expect_null(sim$propensity)
  expect_equal(dim(simulate_exposure(scns$latent, 10, seed = 1,
                                     keep_propensity = TRUE)$propensity), c(10, 2))
})

test_that("N is validated", {
  expect_error(simulate_exposure(scns$independent, 0), "positive integer")
  expect_error(simulate_exposure(scns$independent, 2.5), "positive integer")
  expect_error(simulate_exposure(scns$independent, c(1, 2)), "positive integer")
  expect_error(simulate_exposure(scns$independent, 10, seed = c(1, 2)), "single")
  expect_error(simulate_exposure(list(), 10))
})
