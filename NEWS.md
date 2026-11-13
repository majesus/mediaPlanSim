# mediaPlanSim 0.2.0

* `calibrate_gaussian()` builds a Gaussian-copula scenario from observed reach
  and duplication, with a fast Hermite expansion of the duplication, projection
  on the nearest correlation matrix, a report of unattainable pairs, optional
  dispersion adjustment (`r2 = "auto"`) and regularisation against overfitting
  (`n_eff`/`kappa`, `shrink`, `rank`).
* `compare_to_observed()` tests the compatibility of a scenario with observed
  respondents by parametric bootstrap.
* `scenario_subset()` extracts plans from a larger universe; Gaussian scenarios
  are no longer limited to 16 vehicles for simulation (the true distribution is
  still computed for at most 16).
* `simulate_exposure(insertion_level = TRUE)` returns one 0/1 column per
  insertion; `export_microdata()` writes CSV or Excel.
* The duplication of Gaussian scenarios is computed with the Hermite expansion
  for `abs(rho) <= 0.95`, which is faster and agrees with the previous numerical
  integration to 1e-6.
* Warnings of `qbeta()` for very small shape parameters (rare vehicles) are
  silenced; the reach is preserved.

# mediaPlanSim 0.1.0

* First release.
* Scenarios: `scenario_independent()`, `scenario_gaussian()` and
  `scenario_latent_class()` define data-generating mechanisms for
  multi-insertion exposure.
* `simulate_exposure()` draws synthetic respondents; `exposure_truth()` and
  `population_moments()` return the true distribution of contacts and the
  inputs a model would receive (reach of one and two insertions, duplication);
  `sample_moments()`, `empirical_distribution()` and `distribution_error()` help
  to score a model.
