#' Population moments of a scenario
#'
#' Returns the quantities that a media planner would give to an
#' exposure-distribution model: the reach of one insertion (`r1`) and of two
#' insertions (`r2`) of each vehicle, and the duplication between vehicles.
#'
#' With propensity `P` of a vehicle, `r1 = E[P]`, `r2 = 2 E[P] - E[P^2]` and the
#' duplication of the first insertions of vehicles `i` and `j` is
#' `E[P_i P_j]`. All three are exact for [scenario_independent()] and
#' [scenario_latent_class()]; for [scenario_gaussian()] the reaches are exact
#' and the duplication is computed by deterministic numerical integration with
#' an error far below `1e-6`.
#'
#' @param scenario An object created by one of the `scenario_*()` functions.
#' @return A list with `r1`, `r2` (named numeric vectors), `dup` (symmetric
#'   matrix whose diagonal is `r1`) and `insertions`.
#' @examples
#' scn <- scenario_gaussian(c(A = 0.2, B = 0.15), c(0.312, 0.23925), c(5, 5),
#'                          rho = 0.5)
#' population_moments(scn)
#' @export
population_moments <- function(scenario) {
  stopifnot(inherits(scenario, "mp_scenario"))
  nm <- scenario$names
  m <- scenario$m
  if (scenario$type == "latent_class") {
    w <- scenario$class_prob
    pr <- scenario$class_propensity
    r1 <- as.numeric(crossprod(w, pr))
    e2 <- as.numeric(crossprod(w, pr^2))
    r2 <- 2 * r1 - e2
    d <- crossprod(pr * sqrt(w))
  } else {
    r1 <- scenario$r1
    r2 <- scenario$r2
    d <- outer(r1, r1)
    if (scenario$type == "gaussian" && m > 1L) {
      for (i in seq_len(m - 1L)) for (j in (i + 1L):m) {
        d[i, j] <- d[j, i] <- .dup_gauss(scenario, i, j)
      }
    }
  }
  diag(d) <- r1
  dimnames(d) <- list(nm, nm)
  names(r1) <- names(r2) <- nm
  list(r1 = r1, r2 = r2, dup = d, insertions = stats::setNames(scenario$insertions, nm))
}

# E[g_i(Z1) g_j(rho Z1 + sqrt(1 - rho^2) Z2)] by the trapezoid rule on a normal
# grid, with monotone interpolation of the propensity as a function of z.
.dup_gauss <- function(scn, i, j) {
  rho <- scn$rho[i, j]
  if (abs(rho) < 1e-14) return(scn$r1[i] * scn$r1[j])
  zg <- seq(-10, 10, by = 0.02)
  gi <- splinefun(zg, .pi_from_z(scn$alpha[i], scn$beta[i], zg), method = "monoH.FC")
  gj <- splinefun(zg, .pi_from_z(scn$alpha[j], scn$beta[j], zg), method = "monoH.FC")
  x <- seq(-8, 8, by = 0.025)
  w <- dnorm(x)
  w <- w / sum(w)
  z2 <- outer(rho * x, sqrt(1 - rho^2) * x, "+")
  z2 <- pmax(pmin(z2, 10), -10)
  sum(w * gi(x) * (matrix(gj(z2), nrow(z2)) %*% w))
}

#' True distribution of contacts of a scenario
#'
#' Computes the exact (or numerically integrated) distribution of the total
#' number of contacts of the plan under a scenario, together with its reach,
#' average opportunities to see, gross rating points and the population moments.
#'
#' The distribution is the expectation, over the propensities of a person, of the
#' convolution of the vehicles' binomial distributions. It is evaluated exactly
#' for [scenario_independent()] (convolution of Beta-Binomial distributions) and
#' for [scenario_latent_class()] (finite mixture). For [scenario_gaussian()] the
#' expectation is computed by antithetic Halton quasi-Monte Carlo integration,
#' which is deterministic: it uses no random numbers and does not alter the
#' random state. Its error decreases with `n_qmc`; compare two values of `n_qmc`
#' to judge it. With one or two vehicles the error is typically below `1e-4` for
#' the default `n_qmc`, and it grows with the number of vehicles.
#'
#' @param scenario An object created by one of the `scenario_*()` functions.
#' @param n_qmc Number of quasi-Monte Carlo points for [scenario_gaussian()]
#'   (default 32768); ignored by the other scenarios.
#' @return An object of class `mp_truth`, a list with `distribution` (data
#'   frame with `k`, `pmf`, `cdf` and `reach_k_plus`), `metrics` (`grp`,
#'   `reach_1plus`, `reach_3plus`, `ots`, `mean_contacts`, `variance`),
#'   `moments` (see [population_moments()]), `method` (`"exact"` or `"qmc"`) and
#'   `n_qmc`.
#' @examples
#' scn <- scenario_latent_class(
#'   class_prob = c(0.5, 0.5),
#'   class_propensity = cbind(A = c(0.4, 0.05), B = c(0.05, 0.4)),
#'   insertions = c(4, 4)
#' )
#' exposure_truth(scn)
#' @export
exposure_truth <- function(scenario, n_qmc = 32768L) {
  stopifnot(inherits(scenario, "mp_scenario"))
  n <- scenario$insertions
  m <- scenario$m
  method <- "exact"
  nq <- NA_integer_
  if (scenario$type == "independent") {
    pmf <- 1
    for (j in seq_len(m)) pmf <- .conv_vec(pmf, .betabinom_pmf(n[j], scenario$alpha[j], scenario$beta[j]))
  } else if (scenario$type == "latent_class") {
    pmf <- 0
    for (cc in seq_along(scenario$class_prob)) {
      pc <- 1
      for (j in seq_len(m)) {
        pc <- .conv_vec(pc, dbinom(0:n[j], n[j], scenario$class_propensity[cc, j]))
      }
      pmf <- pmf + scenario$class_prob[cc] * pc
    }
  } else {
    if (length(n_qmc) != 1L || !is.finite(n_qmc) || n_qmc < 256) {
      stop("'n_qmc' must be a single number of at least 256.", call. = FALSE)
    }
    method <- "qmc"
    nq <- as.integer(n_qmc)
    pmf <- .gauss_pmf(scenario, nq)
  }
  dist <- .dist_frame(pmf)
  structure(list(distribution = dist, metrics = .dist_metrics(dist),
                 moments = population_moments(scenario), method = method,
                 n_qmc = nq, scenario = scenario),
            class = "mp_truth")
}

.betabinom_pmf <- function(n, a, b) {
  k <- 0:n
  exp(lchoose(n, k) + lbeta(k + a, n - k + b) - lbeta(a, b))
}

.gauss_pmf <- function(scn, n_qmc) {
  m <- scn$m
  n <- scn$insertions
  half <- ceiling(n_qmc / 2)
  z0 <- qnorm(.halton(half, m))
  z <- rbind(z0, -z0)
  n_sim <- nrow(z)
  if (m > 1L) {
    s_hat <- crossprod(z) / (n_sim - 1L)
    z <- z %*% backsolve(chol(s_hat), diag(m))
  }
  z <- z %*% t(scn$chol_factor)
  acc <- NULL
  for (j in seq_len(m)) {
    pj <- .pi_from_z(scn$alpha[j], scn$beta[j], z[, j])
    pm <- matrix(dbinom(rep(0:n[j], each = n_sim), n[j], rep(pj, n[j] + 1L)),
                 n_sim, n[j] + 1L)
    acc <- if (is.null(acc)) pm else .row_convolve(acc, pm)
  }
  colMeans(acc)
}

#' @export
print.mp_truth <- function(x, digits = 4, max_rows = 12, ...) {
  mt <- x$metrics
  cat("True exposure distribution (mediaPlanSim): ", x$scenario$type,
      " scenario, ", x$method, " evaluation\n", sep = "")
  cat(sprintf("GRPs: %.*f | Reach 1+: %.*f | Reach 3+: %s | OTS: %.*f\n\n",
              2, mt$grp, digits, mt$reach_1plus,
              if (is.na(mt$reach_3plus)) "NA" else sprintf("%.*f", digits, mt$reach_3plus),
              2, mt$ots))
  print(round(x$distribution[seq_len(min(nrow(x$distribution), max_rows)), ], digits),
        row.names = FALSE)
  if (nrow(x$distribution) > max_rows) {
    cat(sprintf("... (%d more rows)\n", nrow(x$distribution) - max_rows))
  }
  invisible(x)
}
