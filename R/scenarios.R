#' Data-generating mechanisms (scenarios)
#'
#' A scenario describes how the propensities of a person to be exposed to each
#' vehicle are generated and how the vehicles depend on each other. Given the
#' propensities, the number of contacts with vehicle `j` is
#' `Binomial(insertions[j], propensity[j])`, and the vehicles are conditionally
#' independent. Within a vehicle, insertions are therefore exchangeable given the
#' propensity.
#'
#' * `scenario_independent()`: propensities are independent across vehicles,
#'   `P_j ~ Beta`, with shape parameters set so that the reach of one insertion
#'   is `r1` and of two insertions `r2` (the contacts of each vehicle are
#'   Beta-Binomial). Duplication is `r1_i * r1_j`.
#' * `scenario_gaussian()`: the same Beta marginals, linked by a Gaussian copula
#'   with correlation `rho` (a number for equicorrelation, or a correlation
#'   matrix). Positive correlation produces duplication above independence, and
#'   negative correlation below it.
#' * `scenario_latent_class()`: the audience is a mixture of segments. Segment
#'   `c` has probability `class_prob[c]` and propensity `class_propensity[c, j]`
#'   for vehicle `j`. Segments that follow different vehicles produce negative
#'   association between them; this is the main way to generate duplication
#'   below independence.
#'
#' For the first two scenarios, `r1` and `r2` are the reach of one and two
#' insertions of each vehicle; `r2` must satisfy `r1 < r2 < 2 * r1 - r1^2`. For
#' the latent-class scenario they are implied by the segments (see
#' [population_moments()]).
#'
#' @param r1,r2 Numeric vectors of reaches (proportions strictly between 0 and
#'   1), one per vehicle. Names give the vehicle names.
#' @param insertions Positive integers, the number of insertions of each
#'   vehicle (at most 400 in total).
#' @param rho A number (equicorrelation) or a symmetric correlation matrix of
#'   the latent Gaussian vector.
#' @param class_prob Numeric vector with the probability of each segment; it is
#'   normalized to sum to one.
#' @param class_propensity Matrix with one row per segment and one column per
#'   vehicle, with probabilities of exposure to one insertion in `[0, 1]`.
#'   Column names give the vehicle names.
#' @return An object of class `mp_scenario`.
#' @seealso [simulate_exposure()], [exposure_truth()], [population_moments()]
#' @examples
#' scn <- scenario_gaussian(r1 = c(A = 0.20, B = 0.15), r2 = c(0.312, 0.23925),
#'                          insertions = c(5, 5), rho = 0.4)
#' scn
#'
#' # Two vehicles with disjoint audiences (negative association)
#' scenario_latent_class(
#'   class_prob = c(0.5, 0.5),
#'   class_propensity = cbind(A = c(0.4, 0.05), B = c(0.05, 0.4)),
#'   insertions = c(4, 4)
#' )
#' @name scenarios
NULL

.new_scenario <- function(type, m, names, n, ...) {
  structure(list(type = type, m = m, names = names, insertions = n, ...),
            class = c(paste0("mp_scenario_", type), "mp_scenario"))
}

#' @rdname scenarios
#' @export
scenario_independent <- function(r1, r2, insertions) {
  mg <- .check_marginals(r1, r2, insertions)
  .new_scenario("independent", mg$m, mg$names, mg$n, r1 = mg$r1, r2 = mg$r2,
                alpha = mg$alpha, beta = mg$beta)
}

#' @rdname scenarios
#' @export
scenario_gaussian <- function(r1, r2, insertions, rho = 0) {
  mg <- .check_marginals(r1, r2, insertions)
  cf <- .corr_factor(rho, mg$m, mg$names)
  .new_scenario("gaussian", mg$m, mg$names, mg$n, r1 = mg$r1, r2 = mg$r2,
                alpha = mg$alpha, beta = mg$beta, rho = cf$R, chol_factor = cf$L)
}

#' Restrict a scenario to some vehicles
#'
#' Extracts a plan from a scenario that describes a larger universe of vehicles,
#' for example one calibrated with [calibrate_gaussian()] on all the titles of a
#' survey. Marginals and dependence of the selected vehicles are kept, and the
#' number of insertions can be changed. The true distribution of a Gaussian
#' scenario is computed for at most 16 vehicles; simulation has no such limit.
#'
#' @param scenario An object created by a `scenario_*()` function.
#' @param vehicles Names or positions of the vehicles to keep.
#' @param insertions Optional new numbers of insertions, one per kept vehicle.
#' @return A `mp_scenario` of the same type.
#' @examples
#' scn <- scenario_gaussian(c(A = 0.2, B = 0.15, C = 0.1),
#'                          c(0.312, 0.23925, 0.17), c(1, 1, 1), rho = 0.3)
#' scenario_subset(scn, c("A", "C"), insertions = c(4, 2))
#' @export
scenario_subset <- function(scenario, vehicles, insertions = NULL) {
  stopifnot(inherits(scenario, "mp_scenario"))
  idx <- if (is.character(vehicles)) match(vehicles, scenario$names) else vehicles
  if (anyNA(idx) || any(idx < 1L | idx > scenario$m) || anyDuplicated(idx)) {
    stop("'vehicles' must be distinct names or positions of the scenario.",
         call. = FALSE)
  }
  n <- if (is.null(insertions)) scenario$insertions[idx] else insertions
  if (length(n) != length(idx)) {
    stop("'insertions' needs one value per kept vehicle.", call. = FALSE)
  }
  nm <- scenario$names[idx]
  switch(
    scenario$type,
    independent = scenario_independent(stats::setNames(scenario$r1[idx], nm),
                                       scenario$r2[idx], n),
    gaussian = scenario_gaussian(stats::setNames(scenario$r1[idx], nm),
                                 scenario$r2[idx], n,
                                 `dimnames<-`(scenario$rho[idx, idx, drop = FALSE],
                                              list(nm, nm))),
    latent_class = scenario_latent_class(
      scenario$class_prob, scenario$class_propensity[, idx, drop = FALSE], n)
  )
}

#' @rdname scenarios
#' @export
scenario_latent_class <- function(class_prob, class_propensity, insertions) {
  pr <- as.matrix(class_propensity)
  if (!is.numeric(class_prob) || !is.numeric(pr) || anyNA(class_prob) ||
      anyNA(pr)) {
    stop("'class_prob' and 'class_propensity' must be numeric without NAs.",
         call. = FALSE)
  }
  if (nrow(pr) != length(class_prob) || nrow(pr) < 1L) {
    stop("'class_propensity' needs one row per element of 'class_prob'.",
         call. = FALSE)
  }
  if (any(class_prob < 0) || sum(class_prob) <= 0) {
    stop("'class_prob' must be non-negative with a positive sum.", call. = FALSE)
  }
  if (any(pr < 0 | pr > 1)) {
    stop("'class_propensity' must lie in [0, 1].", call. = FALSE)
  }
  m <- ncol(pr)
  .check_insertions(insertions, m)
  nm <- colnames(pr)
  if (is.null(nm)) nm <- paste0("V", seq_len(m))
  if (anyDuplicated(nm)) stop("Vehicle names must be unique.", call. = FALSE)
  w <- as.numeric(class_prob) / sum(class_prob)
  dimnames(pr) <- list(NULL, nm)
  .new_scenario("latent_class", m, nm, as.integer(round(insertions)),
                class_prob = w, class_propensity = pr)
}

#' @export
print.mp_scenario <- function(x, ...) {
  cat("Exposure scenario (mediaPlanSim): ", x$type, "\n", sep = "")
  tab <- data.frame(vehicle = x$names, insertions = x$insertions)
  pm <- population_moments(x)
  tab$r1 <- pm$r1
  tab$r2 <- pm$r2
  print(tab, row.names = FALSE, digits = 4)
  if (x$type == "gaussian") {
    cat("Latent correlation:\n")
    print(round(x$rho, 3))
  }
  if (x$type == "latent_class") {
    cat(sprintf("Segments: %d (probabilities %s)\n", length(x$class_prob),
                paste(format(x$class_prob, digits = 3), collapse = ", ")))
  }
  invisible(x)
}
