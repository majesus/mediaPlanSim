# Internal helpers (not exported).

# Beta parameters of a vehicle's exposure propensity P from the reach of one
# (p = r1) and of two insertions (q = r2): E[P] = p and E[P^2] = 2 p - q.
.beta_params <- function(p, q, name) {
  if (!(q > p && q < 2 * p - p^2)) {
    stop(sprintf(paste0("Vehicle '%s': r2 must lie strictly between r1 and ",
                        "2*r1 - r1^2 (in proportions: (%.6g, %.6g)); got %.6g."),
                 name, p, 2 * p - p^2, q), call. = FALSE)
  }
  r11 <- 2 * p - q
  theta <- (r11 - p^2) / (p * (1 - p))
  conc <- (1 - theta) / theta
  c(alpha = p * conc, beta = (1 - p) * conc)
}

.check_marginals <- function(r1, r2, insertions) {
  if (!is.numeric(r1) || !is.numeric(r2) || !is.numeric(insertions)) {
    stop("'r1', 'r2' and 'insertions' must be numeric.", call. = FALSE)
  }
  m <- length(r1)
  if (m < 1L || length(r2) != m || length(insertions) != m) {
    stop("'r1', 'r2' and 'insertions' must be non-empty and of equal length.",
         call. = FALSE)
  }
  if (anyNA(c(r1, r2, insertions)) || !all(is.finite(c(r1, r2, insertions)))) {
    stop("'r1', 'r2' and 'insertions' must be finite.", call. = FALSE)
  }
  if (any(r1 <= 0 | r1 >= 1) || any(r2 <= 0 | r2 >= 1)) {
    stop("'r1' and 'r2' are proportions and must lie strictly between 0 and 1.",
         call. = FALSE)
  }
  .check_insertions(insertions, m)
  nm <- names(r1)
  if (is.null(nm)) nm <- paste0("V", seq_len(m))
  if (anyDuplicated(nm)) stop("Vehicle names must be unique.", call. = FALSE)
  pars <- vapply(seq_len(m), function(i) .beta_params(r1[i], r2[i], nm[i]),
                 c(alpha = 0, beta = 0))
  list(m = m, names = nm, r1 = as.numeric(r1), r2 = as.numeric(r2),
       n = as.integer(round(insertions)), alpha = pars["alpha", ],
       beta = pars["beta", ])
}

.check_insertions <- function(insertions, m) {
  if (length(insertions) != m || !is.numeric(insertions) ||
      anyNA(insertions) || any(insertions < 1) ||
      any(abs(insertions - round(insertions)) > 1e-8)) {
    stop("'insertions' must be positive integers, one per vehicle.",
         call. = FALSE)
  }
  if (sum(insertions) > 400) {
    stop("The total number of insertions is limited to 400.", call. = FALSE)
  }
  invisible(TRUE)
}

# Correlation matrix from a scalar (equicorrelation) or a full matrix, with a
# factor L such that R = L %*% t(L). Positive semi-definite matrices are allowed.
.corr_factor <- function(rho, m, nm) {
  if (m == 1L) {
    r <- matrix(1, 1, 1)
  } else if (is.matrix(rho)) {
    r <- rho
    if (!is.numeric(r) || any(dim(r) != m) || anyNA(r)) {
      stop("'rho' must be a numeric m x m matrix without missing values.",
           call. = FALSE)
    }
    if (max(abs(r - t(r))) > 1e-10) stop("'rho' must be symmetric.", call. = FALSE)
    if (max(abs(diag(r) - 1)) > 1e-10) {
      stop("'rho' must have a unit diagonal.", call. = FALSE)
    }
  } else {
    if (!is.numeric(rho) || length(rho) != 1L || !is.finite(rho)) {
      stop("'rho' must be a single number or a correlation matrix.",
           call. = FALSE)
    }
    r <- matrix(rho, m, m)
    diag(r) <- 1
  }
  ev <- eigen(r, symmetric = TRUE)
  if (min(ev$values) < -1e-8) {
    stop("'rho' is not a valid (positive semi-definite) correlation matrix.",
         call. = FALSE)
  }
  l <- ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), m)
  dimnames(r) <- list(nm, nm)
  list(R = r, L = l)
}

# Propensity as a function of a standard normal latent value (Gaussian copula
# with Beta marginals). The upper-tail form keeps precision in the tails.
.pi_from_z <- function(alpha, beta, z) {
  qbeta(pnorm(z, lower.tail = FALSE), alpha, beta, lower.tail = FALSE)
}

# Run an expression with a given seed and restore the caller's RNG state.
.with_seed <- function(seed, expr) {
  if (is.null(seed)) return(expr)
  if (length(seed) != 1L || !is.finite(seed)) {
    stop("'seed' must be a single finite number or NULL.", call. = FALSE)
  }
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit({
    if (had) assign(".Random.seed", old, envir = globalenv())
    else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  }, add = TRUE)
  set.seed(seed)
  expr
}

# Deterministic shifted Halton points in (0, 1)^dim (no use of the RNG).
.radical_inverse <- function(i, base) {
  out <- numeric(length(i))
  f <- 1 / base
  k <- i
  while (any(k > 0)) {
    out <- out + f * (k %% base)
    k <- k %/% base
    f <- f / base
  }
  out
}

.halton <- function(n, dim) {
  primes <- c(2L, 3L, 5L, 7L, 11L, 13L, 17L, 19L, 23L, 29L, 31L, 37L, 41L,
              43L, 47L, 53L)
  if (dim > length(primes)) {
    stop("At most 16 vehicles are supported with a Gaussian copula.",
         call. = FALSE)
  }
  idx <- seq_len(n)
  shift <- (seq_len(dim) * 0.6180339887498949) %% 1
  u <- vapply(seq_len(dim), function(j) {
    (.radical_inverse(idx, primes[j]) + shift[j]) %% 1
  }, numeric(n))
  u <- matrix(u, nrow = n, ncol = dim)
  pmin(pmax(u, 1e-12), 1 - 1e-12)
}

# Convolution of two probability vectors (support 0, 1, ...).
.conv_vec <- function(a, b) {
  out <- numeric(length(a) + length(b) - 1L)
  lb <- length(b)
  for (i in seq_along(a)) out[i:(i + lb - 1L)] <- out[i:(i + lb - 1L)] + a[i] * b
  out
}

# Row-wise convolution of two matrices whose rows are pmfs.
.row_convolve <- function(a, b) {
  kc <- ncol(a) + ncol(b) - 1L
  out <- matrix(0, nrow(a), kc)
  for (j in seq_len(ncol(a))) {
    cols <- j:(j + ncol(b) - 1L)
    out[, cols] <- out[, cols] + a[, j] * b
  }
  out
}

# Distribution data frame and metrics from a pmf on 0, ..., K.
.dist_frame <- function(pmf) {
  pmf <- pmax(pmf, 0)
  pmf <- pmf / sum(pmf)
  cdf <- pmin(cumsum(pmf), 1)
  data.frame(k = seq_along(pmf) - 1L, pmf = pmf, cdf = cdf,
             reach_k_plus = c(1, 1 - cdf[-length(cdf)]))
}

.dist_metrics <- function(dist) {
  k <- dist$k
  mean_c <- sum(k * dist$pmf)
  reach1 <- 1 - dist$pmf[1]
  list(grp = 100 * mean_c, reach_1plus = reach1,
       reach_3plus = if (nrow(dist) >= 4L) 1 - sum(dist$pmf[1:3]) else NA_real_,
       ots = if (reach1 > 0) mean_c / reach1 else NA_real_,
       mean_contacts = mean_c, variance = sum((k - mean_c)^2 * dist$pmf))
}
