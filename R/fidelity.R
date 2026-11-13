#' Compare a scenario with observed respondents
#'
#' Checks whether the distribution of the total number of contacts observed in a
#' survey is compatible with the one implied by a scenario. The scenario's
#' distribution is computed with [exposure_truth()]. The observed distribution
#' is compared with it by total variation distance, and the distance is judged
#' against a parametric bootstrap: samples of effective size `n_eff` are drawn
#' from the scenario, and their distance to the scenario's distribution is the
#' distance expected from sampling alone if the scenario were the truth. The
#' reported p-value is the share of bootstrap distances at least as large as the
#' observed one.
#'
#' Read the result with these limits in mind. The bootstrap does not
#' re-calibrate the scenario in each replicate, so when the scenario was
#' calibrated on the same respondents the observed distance is optimistic and the
#' test is conservative; calibrate on one part of the sample and compare on the
#' other. Clustering and weighting are summarized in `n_eff` only. A scenario
#' that is not rejected is compatible with the observed distribution, which does
#' not mean that it is the mechanism that generated the data.
#'
#' @param scenario An object created by a `scenario_*()` function or by
#'   [calibrate_gaussian()].
#' @param contacts Observed contacts: a numeric vector with the total number of
#'   contacts of each person, or a matrix with one column per vehicle, in the
#'   same order as the scenario.
#' @param weights Optional non-negative person weights.
#' @param n_eff Effective sample size. Defaults to Kish's
#'   `sum(w)^2 / sum(w^2)`.
#' @param B Number of bootstrap replicates (default 200).
#' @param n_qmc Nodes of the quasi-Monte Carlo integration of the truth.
#' @param seed Optional seed; the global random state is restored.
#' @return A list with `tv_observed`, `tv_bootstrap_mean`, `tv_bootstrap_p95`,
#'   `p_value`, `n_eff`, `observed` and `scenario_distribution` (both
#'   probability vectors on `0, 1, ...`), and the differences in reach
#'   (`reach_difference`) and mean contacts (`mean_difference`), scenario minus
#'   observed.
#' @examples
#' scn <- scenario_gaussian(c(A = 0.2, B = 0.15), c(0.312, 0.23925), c(1, 1), 0.4)
#' obs <- as.data.frame(simulate_exposure(scn, 5000, seed = 1))
#' compare_to_observed(scn, obs[, c("A", "B")], B = 50, seed = 2)$p_value
#' @export
compare_to_observed <- function(scenario, contacts, weights = NULL, n_eff = NULL,
                                B = 200L, n_qmc = 32768L, seed = NULL) {
  stopifnot(inherits(scenario, "mp_scenario"))
  tot <- if (is.matrix(contacts) || is.data.frame(contacts)) {
    rowSums(as.matrix(contacts))
  } else as.numeric(contacts)
  if (anyNA(tot) || any(tot < 0) || any(tot != round(tot))) {
    stop("'contacts' must be non-negative integers without missing values.",
         call. = FALSE)
  }
  kmax <- sum(scenario$insertions)
  if (any(tot > kmax)) {
    stop("'contacts' exceeds the total number of insertions of the scenario.",
         call. = FALSE)
  }
  w <- if (is.null(weights)) rep(1, length(tot)) else as.numeric(weights)
  if (length(w) != length(tot) || anyNA(w) || any(w < 0) || sum(w) <= 0) {
    stop("'weights' must be non-negative, one per person, with positive sum.",
         call. = FALSE)
  }
  obs <- as.numeric(tapply(w, factor(tot, 0:kmax), sum))
  obs[is.na(obs)] <- 0
  obs <- obs / sum(obs)
  if (is.null(n_eff)) n_eff <- sum(w)^2 / sum(w^2)
  th <- exposure_truth(scenario, n_qmc = n_qmc)$distribution$pmf
  tv <- function(a, b) 0.5 * sum(abs(a - b))
  nulo <- .with_seed(seed, vapply(seq_len(B), function(b) {
    sm <- simulate_exposure(scenario, N = max(2L, as.integer(round(n_eff))))
    tv(empirical_distribution(sm)$distribution$pmf, th)
  }, 0))
  t_obs <- tv(obs, th)
  k <- seq_along(th) - 1L
  list(tv_observed = t_obs, tv_bootstrap_mean = mean(nulo),
       tv_bootstrap_p95 = unname(stats::quantile(nulo, 0.95)),
       p_value = (1 + sum(nulo >= t_obs)) / (B + 1), n_eff = n_eff,
       observed = obs, scenario_distribution = th,
       reach_difference = (1 - th[1]) - (1 - obs[1]),
       mean_difference = sum(k * th) - sum(k * obs))
}

#' Export simulated respondents
#'
#' Writes the respondents of a simulation to a CSV or Excel file. CSV can be
#' read from any language, for example with `pandas.read_csv()`. Excel needs the
#' suggested package `writexl` and is limited to 1,048,575 rows.
#'
#' @param sim An object created by [simulate_exposure()].
#' @param file Path of the file. The format is taken from its extension
#'   (`.csv` or `.xlsx`) unless `format` is given.
#' @param format `"csv"` or `"xlsx"`.
#' @param level `"vehicle"` (contacts per vehicle) or `"insertion"` (one 0/1
#'   column per insertion; requires `simulate_exposure(insertion_level = TRUE)`).
#' @param propensity If `TRUE`, the latent propensities are appended (they must
#'   have been kept with `keep_propensity = TRUE`).
#' @return The path of the file, invisibly.
#' @examples
#' sim <- simulate_exposure(scenario_independent(c(A = 0.2), c(0.312), 5), 100,
#'                          seed = 1)
#' f <- tempfile(fileext = ".csv")
#' export_microdata(sim, f)
#' head(read.csv(f), 3)
#' @export
export_microdata <- function(sim, file, format = NULL,
                             level = c("vehicle", "insertion"),
                             propensity = FALSE) {
  stopifnot(inherits(sim, "mp_sim"))
  level <- match.arg(level)
  if (is.null(format)) {
    format <- tolower(tools::file_ext(file))
  }
  if (!format %in% c("csv", "xlsx")) {
    stop("The format must be 'csv' or 'xlsx'.", call. = FALSE)
  }
  df <- as.data.frame(sim, level = level)
  if (propensity) {
    if (is.null(sim$propensity)) {
      stop("The propensities were not kept; use keep_propensity = TRUE.",
           call. = FALSE)
    }
    pp <- as.data.frame(sim$propensity)
    names(pp) <- paste0("propensity_", names(pp))
    df <- cbind(df, pp)
  }
  if (format == "csv") {
    utils::write.csv(df, file, row.names = FALSE)
  } else {
    if (nrow(df) > 1048575L) {
      stop("Excel cannot hold more than 1,048,575 rows; use CSV.", call. = FALSE)
    }
    if (!requireNamespace("writexl", quietly = TRUE)) {
      stop("Package 'writexl' is needed for xlsx; install it or use CSV.",
           call. = FALSE)
    }
    writexl::write_xlsx(df, file)
  }
  invisible(file)
}
