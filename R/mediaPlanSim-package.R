#' mediaPlanSim: Synthetic Multi-Insertion Exposure Data with Known Truth
#'
#' Survey data on media consumption rarely record exposure to several
#' insertions of the same vehicle, so exposure-distribution models cannot be
#' validated against an observed distribution of contacts. This package
#' generates synthetic respondents from explicit data-generating mechanisms
#' (see [scenario_independent()], [scenario_gaussian()] and
#' [scenario_latent_class()]), simulates their contacts with
#' [simulate_exposure()], and returns the true distribution of contacts and the
#' population moments of the mechanism with [exposure_truth()] and
#' [population_moments()]. The inputs of a model (reach of one and two
#' insertions, duplication) can be taken from the truth or estimated from the
#' simulated sample with [sample_moments()], and the model output can be scored
#' with [distribution_error()].
#'
#' **What the package is not.** Data simulated from a mechanism inherit that
#' mechanism's assumptions. A model that reproduces the truth of a mechanism has
#' not thereby been shown to describe real audiences. The package is a tool to
#' study how models behave under controlled and documented conditions.
#'
#' @keywords internal
#' @importFrom stats dbinom dnorm pnorm qbeta qnorm rbinom rnorm runif splinefun
"_PACKAGE"
