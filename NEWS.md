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
