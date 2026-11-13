# Hermite (Mehler) expansion of the duplication of two Gaussian-linked
# propensities: E[g_i(Z1) g_j(rho Z1 + sqrt(1 - rho^2) Z2)] = sum_k rho^k a_ik a_jk,
# where a_ik = E[g_i(Z) He_k(Z)] / sqrt(k!) and He_k are the probabilists'
# Hermite polynomials. The coefficients depend on one vehicle only, so a
# calibration over many pairs reuses them.
.hermite_coef <- function(alpha, beta, K = 400L) {
  x <- seq(-9, 9, by = 0.005)
  w <- dnorm(x)
  w <- w / sum(w)
  g <- .pi_from_z(alpha, beta, x)
  h0 <- rep(1, length(x))
  h1 <- x
  a <- numeric(K + 1L)
  a[1] <- sum(w * g)
  a[2] <- sum(w * g * h1)
  for (k in 2:K) {
    h2 <- (x * h1 - sqrt(k - 1) * h0) / sqrt(k)   # normalized He_k / sqrt(k!)
    a[k + 1L] <- sum(w * g * h2)
    h0 <- h1
    h1 <- h2
  }
  a
}

# Duplication from the coefficients of two vehicles; valid and fast for
# |rho| <= 0.95 (the series converges like |rho|^K).
.dup_series <- function(ai, aj, rho) {
  sum(ai * aj * rho^(0:(length(ai) - 1L)))
}

.rho_bound <- 0.95
.rho_series_max <- 0.995

# Number of Hermite terms needed for a given correlation.
.hermite_terms <- function(rho) {
  rho <- min(abs(rho), .rho_series_max)
  max(400L, min(4000L, as.integer(ceiling(log(1e-9) / log(max(rho, 0.05))))))
}

# Duplication limits of two vehicles over |rho| <= .rho_bound.
.dup_limits <- function(ai, aj) {
  c(lo = .dup_series(ai, aj, -.rho_bound), hi = .dup_series(ai, aj, .rho_bound))
}

# Solve for the correlation that reproduces a target duplication.
.solve_rho <- function(ai, aj, target) {
  lim <- .dup_limits(ai, aj)
  if (target <= lim[["lo"]]) return(c(rho = -.rho_bound, clipped = 1))
  if (target >= lim[["hi"]]) return(c(rho = .rho_bound, clipped = 1))
  f <- function(r) .dup_series(ai, aj, r) - target
  c(rho = stats::uniroot(f, c(-.rho_bound, .rho_bound), tol = 1e-12)$root,
    clipped = 0)
}

# Nearest correlation matrix by alternating projections with Dykstra's
# correction (Higham, 2002).
.nearest_corr <- function(R, tol = 1e-9, maxit = 200L) {
  Y <- R
  dS <- matrix(0, nrow(R), ncol(R))
  for (it in seq_len(maxit)) {
    Rk <- Y - dS
    ev <- eigen(Rk, symmetric = TRUE)
    X <- ev$vectors %*% (pmax(ev$values, 1e-8) * t(ev$vectors))
    dS <- X - Rk
    Yn <- X
    diag(Yn) <- 1
    if (max(abs(Yn - Y)) < tol) { Y <- Yn; break }
    Y <- Yn
  }
  Y <- (Y + t(Y)) / 2
  diag(Y) <- 1
  Y
}

# Low-rank factor correlation: R ~ L L' + diag(1 - rowSums(L^2)).
.factor_corr <- function(R, k, maxit = 100L) {
  m <- nrow(R)
  h <- rep(0.5, m)
  for (it in seq_len(maxit)) {
    Rr <- R
    diag(Rr) <- h
    ev <- eigen(Rr, symmetric = TRUE)
    L <- ev$vectors[, seq_len(k), drop = FALSE] %*%
      diag(sqrt(pmax(ev$values[seq_len(k)], 0)), k)
    hn <- pmin(rowSums(L^2), 0.999)
    if (max(abs(hn - h)) < 1e-8) { h <- hn; break }
    h <- hn
  }
  Rl <- L %*% t(L)
  diag(Rl) <- 1
  Rl
}

#' Calibrate a Gaussian-copula scenario to observed moments
#'
#' Builds a [scenario_gaussian()] whose vehicles reproduce the reach of one
#' insertion (`r1`) and the duplication matrix (`dup`) of an observed audience
#' study, so that synthetic respondents can be generated with the dependence
#' of real data. It is the step that is missing when only aggregate figures are
#' available: for each pair of vehicles the latent correlation that reproduces
#' the observed duplication is found with a fast Hermite expansion, the
#' resulting matrix is projected on the nearest valid correlation matrix when
#' necessary, and every adjustment is reported in `attr(, "calibration")`.
#'
#' A Gaussian copula with Beta marginals cannot reproduce every duplication.
#' Pairs whose duplication lies outside the attainable range for
#' `abs(rho) <= 0.95` are clipped to the limit and listed in
#' `attr(, "calibration")$flagged`, with the absolute error left. A duplication
#' of exactly zero, which a sample of rare vehicles can produce, is always
#' flagged because positive propensities cannot give it. The attainable range widens when the
#' propensities are more dispersed, which corresponds to a smaller `r2`.
#' With `r2 = "auto"` the dispersion of the vehicles involved in a flagged pair
#' is increased stepwise (down to 5 per cent of the admissible range above
#' `r1`) until all pairs are attainable or the limit is reached. This changes
#' `r2`, which is then a modelling assumption and not an observation; the
#' values used are returned in `attr(, "calibration")$r2_used`. If `r2` is
#' `NULL`, the midpoint of the admissible range `(r1, 2 r1 - r1^2)` is used.
#'
#' **Overfitting.** A free correlation for every pair has `m (m - 1) / 2`
#' parameters and reproduces the sampling noise of the duplication matrix.
#' When the effective sample size `n_eff` of the study is given, the
#' correlation of each pair is multiplied by `n_eff * dup / (n_eff * dup + kappa)`,
#' that is, by the share of its expected number of joint exposures in
#' `kappa` more: pairs with few joint exposures, whose duplication is mostly
#' sampling noise (a duplication of zero among rare vehicles, for instance), are
#' pulled towards independence. Besides that, two regularisers are
#' available: `rank` replaces the matrix by a factor
#' structure with that many factors, and `shrink` in `[0, 1]` moves the
#' correlations towards their average. Choose them by calibrating on part of the
#' respondents and checking the fidelity on the rest with
#' [compare_to_observed()]. Duplications estimated from a sample are noisy for
#' rare vehicles, so a calibration to them inherits that noise. The regularisers
#' are off by default because they are not free: in a held-out comparison with
#' one day of newspaper readership (22 titles, strong clustering of titles by
#' publisher), the unregularised calibration was the closest to the held-out
#' respondents and each regulariser moved the scenario towards independence and
#' away from them. With few respondents or many rare vehicles the balance can be
#' the opposite, which is why the choice should be checked on held-out data.
#'
#' @param r1 Reach of one insertion of each vehicle (named numeric vector).
#' @param r2 `NULL`, `"auto"` or the reach of two insertions of each vehicle.
#' @param insertions Number of insertions of each vehicle of the scenario to
#'   generate.
#' @param dup Symmetric matrix of duplication of the first insertions (the
#'   diagonal is ignored and replaced by `r1`).
#' @param rank Optional integer, number of factors of the correlation matrix.
#' @param shrink Number in `[0, 1]`, weight of the average correlation.
#' @param n_eff Optional effective sample size of the study that gave `dup`; it
#'   activates the shrinkage of rare pairs described above.
#' @param kappa Non-negative number, joint exposures that a pair needs to keep
#'   half of its correlation (default 10). Use [compare_to_observed()] on held-out
#'   respondents to choose it.
#' @return A `mp_scenario` of type `gaussian` with an attribute `calibration`:
#'   a list with `target_dup`, `achieved_dup`, `rho_raw` (pair-wise solution),
#'   `rho` (matrix used), `flagged` (data frame of clipped pairs), `r2_used`,
#'   `psd_adjusted` (logical) and `max_abs_error` and `mean_rel_error` of the
#'   achieved duplication.
#' @seealso [compare_to_observed()], [scenario_gaussian()]
#' @examples
#' r1 <- c(A = 0.20, B = 0.15, C = 0.10)
#' dup <- matrix(c(0.20, 0.06, 0.03, 0.06, 0.15, 0.02, 0.03, 0.02, 0.10), 3)
#' scn <- calibrate_gaussian(r1, r2 = NULL, insertions = c(1, 1, 1), dup = dup)
#' attr(scn, "calibration")$max_abs_error
#' @export
calibrate_gaussian <- function(r1, r2 = NULL, insertions, dup, rank = NULL,
                               shrink = 0, n_eff = NULL, kappa = 10) {
  m <- length(r1)
  if (is.null(names(r1))) names(r1) <- paste0("V", seq_len(m))
  nm <- names(r1)
  if (!is.matrix(dup) || any(dim(dup) != m) || anyNA(dup) ||
      max(abs(dup - t(dup))) > 1e-8) {
    stop("'dup' must be a symmetric m x m matrix without missing values.",
         call. = FALSE)
  }
  if (!is.numeric(shrink) || length(shrink) != 1L || shrink < 0 || shrink > 1) {
    stop("'shrink' must be a number in [0, 1].", call. = FALSE)
  }
  if (!is.null(n_eff) && (!is.numeric(n_eff) || length(n_eff) != 1L || n_eff <= 0)) {
    stop("'n_eff' must be a single positive number.", call. = FALSE)
  }
  if (!is.numeric(kappa) || length(kappa) != 1L || kappa < 0) {
    stop("'kappa' must be a single non-negative number.", call. = FALSE)
  }
  if (!is.null(rank) && (length(rank) != 1L || rank < 1 || rank > m ||
                         rank != round(rank))) {
    stop("'rank' must be an integer between 1 and the number of vehicles.",
         call. = FALSE)
  }
  auto <- identical(r2, "auto")
  frac <- rep(0.5, m)
  mid <- function(f) r1 + f * (r1 - r1^2)
  if (is.null(r2) || auto) r2v <- mid(frac) else r2v <- as.numeric(r2)
  names(r2v) <- nm
  d <- dup
  diag(d) <- r1
  dimnames(d) <- list(nm, nm)
  coefs <- function(r2v) lapply(seq_len(m), function(i) {
    p <- .beta_params(unname(r1[i]), unname(r2v[i]), nm[i])
    .hermite_coef(p[["alpha"]], p[["beta"]])
  })
  solve_all <- function(co) {
    rho <- diag(m)
    cl <- matrix(0, m, m)
    for (i in seq_len(m - 1L)) for (j in (i + 1L):m) {
      s <- .solve_rho(co[[i]], co[[j]], d[i, j])
      rho[i, j] <- rho[j, i] <- s[["rho"]]
      cl[i, j] <- cl[j, i] <- s[["clipped"]]
    }
    list(rho = rho, clipped = cl)
  }
  co <- coefs(r2v)
  sol <- solve_all(co)
  if (auto) {
    steps <- c(0.5, 0.35, 0.25, 0.15, 0.1, 0.05)
    lvl <- rep(1L, m)
    repeat {
      bad <- which(rowSums(sol$clipped) > 0)
      bad <- bad[lvl[bad] < length(steps)]
      if (!length(bad)) break
      lvl[bad] <- lvl[bad] + 1L
      frac <- steps[lvl]
      r2v <- stats::setNames(mid(frac), nm)
      co <- coefs(r2v)
      sol <- solve_all(co)
    }
  }
  rho_raw <- sol$rho
  dimnames(rho_raw) <- list(nm, nm)
  R <- rho_raw
  if (!is.null(n_eff) && m > 1L) {
    lam <- (n_eff * d) / (n_eff * d + kappa)
    diag(lam) <- 1
    R <- R * lam
  }
  if (shrink > 0 && m > 1L) {
    avg <- mean(R[upper.tri(R)])
    tgt <- matrix(avg, m, m)
    diag(tgt) <- 1
    R <- (1 - shrink) * R + shrink * tgt
  }
  if (!is.null(rank) && rank < m) R <- .factor_corr(R, as.integer(rank))
  psd <- FALSE
  if (m > 1L && min(eigen(R, symmetric = TRUE, only.values = TRUE)$values) < 1e-8) {
    R <- .nearest_corr(R)
    psd <- TRUE
  }
  dimnames(R) <- list(nm, nm)
  scn <- scenario_gaussian(r1, r2v, insertions, R)
  ach <- population_moments(scn)$dup
  off <- upper.tri(d)
  fl <- which(sol$clipped == 1 & upper.tri(sol$clipped), arr.ind = TRUE)
  flagged <- data.frame(
    vehicle_i = nm[fl[, 1]], vehicle_j = nm[fl[, 2]],
    target = d[fl], achieved = ach[fl],
    stringsAsFactors = FALSE)
  flagged$abs_error <- abs(flagged$achieved - flagged$target)
  attr(scn, "calibration") <- list(
    target_dup = d, achieved_dup = ach, rho_raw = rho_raw, rho = R,
    flagged = flagged, r2_used = r2v, psd_adjusted = psd,
    max_abs_error = if (m > 1L) max(abs(ach - d)[off]) else 0,
    mean_rel_error = if (m > 1L && any(d[off] > 0)) {
      mean((abs(ach - d) / d)[off & d > 0])
    } else 0)
  scn
}
