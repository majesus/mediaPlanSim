r1 <- c(A = 0.20, B = 0.15)
r2 <- c(A = 0.312, B = 0.23925)

test_that("one vehicle reproduces r1 and r2 exactly", {
  for (n in 1:2) {
    tr <- exposure_truth(scenario_independent(c(V = 0.2), c(V = 0.312), n))
    expect_equal(tr$metrics$reach_1plus, c(0.2, 0.312)[n], tolerance = 1e-12)
  }
  tg <- exposure_truth(scenario_gaussian(c(V = 0.2), c(V = 0.312), 2))
  expect_equal(tg$metrics$reach_1plus, 0.312, tolerance = 1e-3)
})

test_that("independent truth equals the product of Beta-Binomial marginals", {
  s <- scenario_independent(r1, r2, c(5, 4))
  tr <- exposure_truth(s)
  a <- s$alpha
  b <- s$beta
  bb <- function(n, a, b) choose(n, 0:n) * beta(0:n + a, n - 0:n + b) / beta(a, b)
  ref <- as.numeric(tapply(outer(bb(5, a[1], b[1]), bb(4, a[2], b[2])),
                           outer(0:5, 0:4, "+"), sum))
  expect_equal(tr$distribution$pmf, ref, tolerance = 1e-12)
  expect_equal(sum(tr$distribution$pmf), 1)
  expect_equal(tr$metrics$mean_contacts, 5 * 0.2 + 4 * 0.15, tolerance = 1e-12)
  expect_equal(tr$metrics$grp, 100 * (5 * 0.2 + 4 * 0.15), tolerance = 1e-10)
  expect_equal(tr$method, "exact")
})

test_that("a single latent class is a plain binomial convolution", {
  s <- scenario_latent_class(1, cbind(A = 0.3, B = 0.3), c(2, 3))
  tr <- exposure_truth(s)
  expect_equal(tr$distribution$pmf, dbinom(0:5, 5, 0.3), tolerance = 1e-12)
})

test_that("latent-class moments follow the mixture formulas", {
  w <- c(0.35, 0.35, 0.30)
  p <- cbind(A = c(.45, .15, .15), B = c(.15, .45, .15))
  pm <- population_moments(scenario_latent_class(w, p, c(4, 4)))
  expect_equal(unname(pm$r1["A"]), sum(w * p[, 1]), tolerance = 1e-14)
  expect_equal(unname(pm$r2["A"]), 2 * sum(w * p[, 1]) - sum(w * p[, 1]^2),
               tolerance = 1e-14)
  expect_equal(unname(pm$dup["A", "B"]), sum(w * p[, 1] * p[, 2]), tolerance = 1e-14)
  expect_lt(pm$dup["A", "B"], prod(pm$r1))        # negative association
})

test_that("Gaussian copula: rho = 0 gives independence, sign orders duplication", {
  s0 <- scenario_gaussian(r1, r2, c(5, 5), rho = 0)
  si <- scenario_independent(r1, r2, c(5, 5))
  expect_lt(distribution_error(exposure_truth(s0), exposure_truth(si))[["tv"]], 2e-3)
  d <- sapply(c(-0.4, 0, 0.4, 0.8), function(r) {
    population_moments(scenario_gaussian(r1, r2, c(5, 5), rho = r))$dup[1, 2]
  })
  expect_true(all(diff(d) > 0))
  expect_equal(d[2], 0.20 * 0.15, tolerance = 1e-12)
})

test_that("Gaussian duplication matches a large simulation", {
  s <- scenario_gaussian(c(A = .2, B = .15, C = .3), c(.312, .23925, .45),
                         c(5, 5, 4), rho = 0.5)
  pm <- population_moments(s)
  sim <- simulate_exposure(s, 400000, seed = 11, keep_propensity = TRUE)
  d_mc <- crossprod(sim$propensity) / sim$N
  expect_lt(max(abs((d_mc - pm$dup)[upper.tri(d_mc)])), 1.5e-3)
})

test_that("quasi-Monte Carlo truth is stable in n_qmc and the RNG is untouched", {
  s <- scenario_gaussian(c(A = .2, B = .15, C = .3), c(.312, .23925, .45),
                         c(5, 5, 4), rho = 0.5)
  set.seed(1)
  before <- runif(1)
  set.seed(1)
  a <- exposure_truth(s, 8192)
  b <- exposure_truth(s, 65536)
  after <- runif(1)
  expect_equal(before, after)
  expect_lt(distribution_error(a, b)[["tv"]], 2e-3)
  expect_equal(a$method, "qmc")
  expect_equal(a$n_qmc, 8192L)
  expect_error(exposure_truth(s, 10), "at least 256")
})

test_that("truth print works and pmf is valid", {
  tr <- exposure_truth(scenario_independent(r1, r2, c(5, 5)))
  expect_output(print(tr), "GRPs")
  expect_true(all(tr$distribution$pmf >= 0))
  expect_equal(tr$distribution$cdf[nrow(tr$distribution)], 1)
  expect_equal(tr$distribution$reach_k_plus[1], 1)
})
