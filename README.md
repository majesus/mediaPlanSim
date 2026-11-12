# mediaPlanSim

Synthetic multi-insertion exposure data **with a known truth** for media plans.

Survey data rarely record exposure to several insertions of the same vehicle, so
models of the distribution of contacts cannot be checked against an observed
distribution. `mediaPlanSim` generates synthetic respondents from explicit
mechanisms and returns the true distribution of contacts of the plan together
with the inputs a planner would give to a model: the reach of one and of two
insertions of each vehicle and the duplication between vehicles.

## Installation

```r
# install.packages("remotes")
remotes::install_github("majesus/mediaPlanSim")
```

## Example

```r
library(mediaPlanSim)

scn <- scenario_gaussian(
  r1 = c(TV = 0.20, Press = 0.15), r2 = c(0.312, 0.23925),
  insertions = c(5, 4), rho = 0.4
)

exposure_truth(scn)               # true distribution of contacts, reach, OTS
population_moments(scn)           # r1, r2 and duplication (inputs of a model)

sim <- simulate_exposure(scn, N = 10000, seed = 1)
sample_moments(sim)               # the same moments, estimated from the sample
distribution_error(empirical_distribution(sim), exposure_truth(scn))
write.csv(as.data.frame(sim), "respondents.csv", row.names = FALSE)  # for Python
```

## Mechanisms

| Function | Dependence between vehicles |
|---|---|
| `scenario_independent()` | none; Beta-Binomial marginals |
| `scenario_gaussian()` | Gaussian copula on exposure propensities, positive or negative |
| `scenario_latent_class()` | audience segments; can produce duplication below independence |

Given the propensities of a person, the contacts with each vehicle are binomial
and vehicles are conditionally independent. The truth is exact for the independent
and latent-class mechanisms and deterministic quasi-Monte Carlo for the Gaussian
copula.

## Limits

* Data simulated from a mechanism inherit its assumptions. A model that reproduces
  the truth of a mechanism has not been shown to describe real audiences.
* The reach of one and two insertions and the duplication do not determine the
  distribution of contacts. The vignette shows two mechanisms with the same
  moments and different distributions; benchmarks should use several mechanisms.
* Insertions of a vehicle are exchangeable given the propensity, and only
  vehicles with a bounded number of insertions are generated.

See the vignette `vignette("benchmarking", package = "mediaPlanSim")`.

## License

MIT.
