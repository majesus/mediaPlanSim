r1 <- c(A = 0.14, B = 0.05, C = 0.02)
r2 <- r1 + 0.4 * (r1 - r1^2)
true_rho <- matrix(c(1, .6, .3, .6, 1, .5, .3, .5, 1), 3)

test_that("the Hermite series reproduces the numerical duplication", {
  p <- .beta_params(0.14, 0.14 + 0.4 * (0.14 - 0.14^2), "A")
  q <- .beta_params(0.05, 0.05 + 0.4 * (0.05 - 0.05^2), "B")
  ai <- .hermite_coef(p[["alpha"]], p[["beta"]])
  aj <- .hermite_coef(q[["alpha"]], q[["beta"]])
  scn <- scenario_gaussian(c(A = 0.14, B = 0.05),
                           c(0.14 + 0.4 * (0.14 - 0.14^2), 0.05 + 0.4 * (0.05 - 0.05^2)),
                           c(1, 1), 0.5)
  zg <- seq(-10, 10, by = 0.02)
  gi <- splinefun(zg, .pi_from_z(p[["alpha"]], p[["beta"]], zg), method = "monoH.FC")
  gj <- splinefun(zg, .pi_from_z(q[["alpha"]], q[["beta"]], zg), method = "monoH.FC")
  x <- seq(-8, 8, by = 0.025)
  w <- dnorm(x) / sum(dnorm(x))
  z2 <- pmax(pmin(outer(0.5 * x, sqrt(1 - 0.25) * x, "+"), 10), -10)
  num <- sum(w * gi(x) * (matrix(gj(z2), nrow(z2)) %*% w))
  expect_equal(.dup_series(ai, aj, 0.5), num, tolerance = 1e-6)
  expect_equal(population_moments(scn)$dup[1, 2], num, tolerance = 1e-6)
})

test_that("calibration recovers a known correlation matrix", {
  scn <- scenario_gaussian(r1, r2, c(1, 1, 1), true_rho)
  d <- population_moments(scn)$dup
  cal <- calibrate_gaussian(r1, r2, c(1, 1, 1), d)
  info <- attr(cal, "calibration")
  expect_s3_class(cal, "mp_scenario_gaussian")
  expect_equal(unname(info$rho), unname(true_rho), tolerance = 1e-6)
  expect_lt(info$max_abs_error, 1e-8)
  expect_equal(nrow(info$flagged), 0L)
  expect_false(info$psd_adjusted)
})

test_that("unattainable duplications are flagged and r2 = 'auto' can fix them", {
  scn <- scenario_gaussian(r1, r2, c(1, 1, 1), true_rho)
  d <- population_moments(scn)$dup
  d[1, 2] <- d[2, 1] <- 0.049
  fixed <- calibrate_gaussian(r1, r2, c(1, 1, 1), d)
  expect_equal(nrow(attr(fixed, "calibration")$flagged), 1L)
  expect_lt(attr(fixed, "calibration")$flagged$achieved, 0.049)
  auto <- calibrate_gaussian(r1, "auto", c(1, 1, 1), d)
  expect_equal(nrow(attr(auto, "calibration")$flagged), 0L)
  expect_lt(attr(auto, "calibration")$max_abs_error, 1e-7)
  mid <- r1 + 0.5 * (r1 - r1^2)
  expect_true(any(attr(auto, "calibration")$r2_used < mid - 1e-9))
  expect_true(all(attr(auto, "calibration")$r2_used > r1))
})

test_that("regularisers keep a valid correlation and change the matrix", {
  scn <- scenario_gaussian(r1, r2, c(1, 1, 1), true_rho)
  d <- population_moments(scn)$dup
  sh <- attr(calibrate_gaussian(r1, r2, c(1, 1, 1), d, shrink = 0.5), "calibration")$rho
  rk <- attr(calibrate_gaussian(r1, r2, c(1, 1, 1), d, rank = 1), "calibration")$rho
  for (R in list(sh, rk)) {
    expect_gte(min(eigen(R, symmetric = TRUE, only.values = TRUE)$values), -1e-8)
    expect_equal(unname(diag(R)), c(1, 1, 1))
  }
  expect_false(isTRUE(all.equal(unname(sh), unname(true_rho))))
  expect_error(calibrate_gaussian(r1, r2, c(1, 1, 1), d, shrink = 2), "shrink")
  expect_error(calibrate_gaussian(r1, r2, c(1, 1, 1), d, rank = 7), "rank")
  expect_error(calibrate_gaussian(r1, r2, c(1, 1, 1), d[1:2, 1:2]), "dup")
})

test_that("an inconsistent duplication matrix is projected and reported", {
  d <- matrix(c(0.14, 0.04, 0.001, 0.04, 0.05, 0.04, 0.001, 0.04, 0.02), 3)
  cal <- calibrate_gaussian(r1, r2, c(1, 1, 1), d)
  expect_true(is.logical(attr(cal, "calibration")$psd_adjusted))
  expect_gte(min(eigen(attr(cal, "calibration")$rho, only.values = TRUE)$values),
             -1e-8)
})

test_that("scenario_subset keeps marginals and dependence", {
  scn <- scenario_gaussian(r1, r2, c(1, 1, 1), true_rho)
  sb <- scenario_subset(scn, c("A", "C"), insertions = c(4, 2))
  expect_equal(sb$names, c("A", "C"))
  expect_equal(sb$insertions, c(4L, 2L))
  expect_equal(unname(sb$rho[1, 2]), 0.3)
  expect_equal(unname(population_moments(sb)$r1), unname(r1[c(1, 3)]))
  lc <- scenario_latent_class(c(.5, .5), cbind(A = c(.4, .05), B = c(.05, .4), C = c(.1, .1)),
                              c(1, 1, 1))
  expect_equal(scenario_subset(lc, 2:3)$names, c("B", "C"))
  expect_error(scenario_subset(scn, "Z"), "vehicles")
  big <- scenario_gaussian(setNames(rep(0.1, 18), paste0("V", 1:18)),
                           rep(0.17, 18), rep(1, 18), 0.2)
  expect_error(exposure_truth(big), "16")
  expect_s3_class(simulate_exposure(big, 50, seed = 1), "mp_sim")
})

test_that("compare_to_observed accepts the true scenario and rejects independence", {
  scn <- scenario_gaussian(r1, r2, c(1, 1, 1), true_rho)
  obs <- simulate_exposure(scn, 20000, seed = 11)
  ok <- compare_to_observed(scn, obs$contacts, B = 60, seed = 1)
  bad <- compare_to_observed(scenario_independent(r1, r2, c(1, 1, 1)),
                             obs$contacts, B = 60, seed = 1)
  expect_gt(ok$p_value, 0.01)
  expect_lt(bad$p_value, 0.05)
  expect_lt(ok$tv_observed, bad$tv_observed)
  expect_equal(sum(ok$observed), 1)
  expect_error(compare_to_observed(scn, c(0, 9)), "exceeds")
  expect_error(compare_to_observed(scn, c(0, 1), weights = c(1, -1)), "weights")
})

test_that("insertion-level simulation is consistent with contacts", {
  scn <- scenario_gaussian(r1, r2, c(3, 2, 1), true_rho)
  sim <- simulate_exposure(scn, 4000, seed = 5, insertion_level = TRUE)
  expect_equal(dim(sim$insertions), c(4000L, 6L))
  expect_true(all(sim$insertions %in% 0:1))
  expect_equal(unname(sim$contacts[, 1]), unname(rowSums(sim$insertions[, 1:3])))
  expect_equal(colnames(sim$insertions), c("A_1", "A_2", "A_3", "B_1", "B_2", "C_1"))
  df <- as.data.frame(sim, level = "insertion")
  expect_equal(names(df), c("id", colnames(sim$insertions), "total"))
  expect_error(as.data.frame(simulate_exposure(scn, 10, seed = 1), level = "insertion"),
               "insertion-level")
  m <- sample_moments(sim)
  expect_equal(unname(m$r1), unname(r1), tolerance = 0.15)
})

test_that("export_microdata writes CSV and refuses wrong requests", {
  scn <- scenario_gaussian(r1, r2, c(3, 2, 1), true_rho)
  sim <- simulate_exposure(scn, 50, seed = 5, keep_propensity = TRUE)
  f <- tempfile(fileext = ".csv")
  export_microdata(sim, f, propensity = TRUE)
  d <- read.csv(f)
  expect_equal(nrow(d), 50L)
  expect_true(all(c("id", "A", "B", "C", "total", "propensity_A") %in% names(d)))
  expect_error(export_microdata(sim, tempfile(fileext = ".txt")), "format")
  expect_error(export_microdata(sim, f, level = "insertion"), "insertion-level")
  expect_error(export_microdata(simulate_exposure(scn, 5, seed = 1), f,
                                propensity = TRUE), "propensities")
  skip_if_not_installed("writexl")
  g <- tempfile(fileext = ".xlsx")
  export_microdata(sim, g)
  expect_true(file.exists(g))
})

test_that("rare vehicles keep their reach despite extreme Beta shapes", {
  for (p in c(0.001, 0.01)) {
    pr <- .beta_params(p, p + 0.1 * (p - p^2), "x")
    a <- .hermite_coef(pr[["alpha"]], pr[["beta"]])
    expect_equal(a[1] / p, 1, tolerance = 1e-4)
  }
  scn <- scenario_gaussian(c(A = 0.001, B = 0.01),
                           c(0.001 + 0.2 * (0.001 - 1e-6), 0.01 + 0.2 * (0.01 - 1e-4)),
                           c(1, 1), 0.5)
  expect_silent(sim <- simulate_exposure(scn, 2000, seed = 1))
})

test_that("n_eff shrinks the correlation of pairs with few joint exposures", {
  scn <- scenario_gaussian(r1, r2, c(1, 1, 1), true_rho)
  d <- population_moments(scn)$dup
  big <- attr(calibrate_gaussian(r1, r2, c(1, 1, 1), d, n_eff = 1e9), "calibration")$rho
  small <- attr(calibrate_gaussian(r1, r2, c(1, 1, 1), d, n_eff = 200, kappa = 10), "calibration")$rho
  expect_equal(unname(big), unname(true_rho), tolerance = 1e-5)
  expect_true(all(abs(small[upper.tri(small)]) < abs(true_rho[upper.tri(true_rho)])))
  d0 <- d
  d0[1, 3] <- d0[3, 1] <- 0
  z <- attr(calibrate_gaussian(r1, r2, c(1, 1, 1), d0, n_eff = 5000), "calibration")
  expect_lt(abs(z$rho[1, 3]), 0.2)
  expect_error(calibrate_gaussian(r1, r2, c(1, 1, 1), d, n_eff = -1), "n_eff")
  expect_error(calibrate_gaussian(r1, r2, c(1, 1, 1), d, kappa = -1), "kappa")
})

test_that("extremely dispersed propensities give finite truths", {
  r1x <- c(A = 0.14, B = 0.005)
  scn <- scenario_gaussian(r1x, r1x + 0.02 * (r1x - r1x^2), c(2, 2), 0.4)
  expect_true(all(is.finite(.pi_from_z(scn$alpha[2], scn$beta[2], seq(-9, 9, by = 0.5)))))
  tr <- exposure_truth(scn)
  expect_true(all(is.finite(tr$distribution$pmf)))
  expect_equal(sum(tr$distribution$pmf), 1, tolerance = 1e-8)
})
