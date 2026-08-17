# StateSpaceSIR

An R package for fitting Bayesian state-space generalized logistic (surplus
production) models by sampling-importance-resampling (SIR). The SIR algorithm
follows McAllister et al. (1994); the model structure follows Zerbini et al.
(2011, 2019) and Romero et al. (2022). Version 1.1 adds four parameterisations
of a demographic Allee effect.

## Installation

```r
devtools::install_github("SW-atlantic-right-whale-assessment/StateSpaceSIR")
```

Needs a C++ toolchain (Rtools on Windows, Xcode command line tools on macOS).

## Population dynamics

Abundance in year *y* is projected as

    N[y+1] = (N[y] + f(N[y]) * h(N[y]) - C[y] * SLR[y]) * exp(e[y])

where *f* is the generalized theta-logistic surplus production function

    f(N) = r_max * N * (1 - (N / K)^z)

*h* is a depensation function (below), *C* is annual catch, *SLR* a
period-specific struck-and-loss correction, and *e* lognormal process error with
variance sigma^2. The shape parameter *z* is solved numerically from `Pmsy`, the
proportion of *K* at which maximum production occurs.

*K* is not given a prior. Abundance is projected "backwards": a prior is placed
on a recent abundance `N_recent` and the trajectory is back-calculated. Priors
are therefore specified for `N_recent`, `r_max`, `Pmsy` (or `z`), `SLR`, the catch
parameter, and `sigma^2`. Catch in year *y* is

    C[y] = C_min[y] + theta * (C_max[y] - C_min[y])

with `theta` estimated. Absolute and relative abundance are fitted with lognormal
likelihoods; catchability is integrated out analytically, including for the
multivariate lognormal case with a full among-year variance–covariance matrix.

## Depensation

`allee_model` selects the function *h* applied to surplus production. `P50`
(written *Pd* in Romero et al., in review) is the depletion level controlling the
strength of the effect.

| `allee_model` | `h(N)` | Source |
|---|---|---|
| 0 | `1` | no depensation |
| 1 | `1 - exp(log(0.5) * N / (P50 * K))` | Hilborn et al. (2014) |
| 2 | `2 / (1 + exp(-N / (P50 * K))) - 1` | logistic |
| 3 | `(N / K) - P50` | Lin & Li (2002) |
| 4 | `(P50*K / (K - P50*K)) * (N / (P50*K) - 1)` | Haider et al. (2017) |

Models 1 and 2 are weak Allee effects: per-capita growth stays positive. Models 3
and 4 are strong: growth turns negative below a threshold.

## Example

Fit a model without depensation, plus its prior predictive check:

```r
library(StateSpaceSIR)

sir <- list()
for (i in 1:2) {                 # i = 1 fit; i = 2 prior predictive
  sir[[i]] <- StateSpaceSIR(
    file_name   = NULL,
    allee_model = 0,
    n_resamples = 20000,
    priors = make_prior_list(
      r_max = make_prior(runif, 0, 0.11),
      N_obs = make_prior(runif, 100, 10000),
      var_N = make_prior(runif, 6.506055e-05, 6.506055e-04),
      z     = make_prior(use = FALSE),
      Pmsy  = make_prior(runif, 0.5, 0.8)),
    catch_multipliers = make_multiplier_list(   # one per catch period
      make_prior(1),
      make_prior(rnorm, 1.60, 0.04),
      make_prior(rnorm, 1.09, 0.04),
      make_prior(1)),
    target.Yr        = 2019,
    num.haplotypes   = 24,
    output.Yrs       = c(2021, 2030),
    abs.abundance    = Abs.Abundance.2005,
    abs.abundance.key = TRUE,
    rel.abundance    = Rel.Abundance,
    rel.abundance.key = TRUE,
    count.data       = Count.Data,
    count.data.key   = FALSE,
    growth.rate.obs  = c(0.074, 0.033, FALSE),
    growth.rate.Yrs  = c(1995, 1996, 1997, 1998),
    catch.data       = Catch.data,
    control          = sir_control(threshold = 1e-5, progress_bar = TRUE),
    realized_prior   = (i == 2))
}

summary_sir(sir[[1]]$resamples_output, object = "Resample_Summary")
plot_trajectory(sir[[1]])
plot_density(SIR = list(sir[[1]]), priors = list(sir[[2]]), inc_reference = FALSE)
plot_ioa(sir[[1]])
summary_table(sir[[1]])
```

To add depensation, set `allee_model` to 1–4 and give `P50` a prior in
`make_prior_list()`.

`Catch.data`, `Rel.Abundance`, `Abs.Abundance.2005` and `Count.Data` ship with
the package (`?Catch.data`).

## Model comparison

`bayes_factor()` and `waic()` take a list of fitted models; `weight_model()`
builds a model-averaged posterior by resampling in proportion to Bayes factors.
Models must share a likelihood to be comparable.

## Applications

* Romero et al. (2022), southern right whales, southwestern Atlantic —
  [code and data](https://github.com/SW-atlantic-right-whale-2021-assessment)
* Romero et al. (in review, *Ecology*), depensation in southern right whales —
  [DepensationRuns](https://github.com/SW-atlantic-right-whale-assessment/DepensationRuns)

## Citation

`citation("StateSpaceSIR")`, or see `CITATION.cff`. The Zenodo DOI is added on
release.

## License

MIT. See `LICENSE.md`.

## References

Haider, H. S., et al. 2017. Incorporating Allee effects into the potential
biological removal level. *Natural Resource Modeling* 30:e12133.

Hilborn, R., D. J. Hively, O. P. Jensen, and T. A. Branch. 2014. The dynamics of
fish populations at low abundance and prospects for rebuilding and recovery.
*ICES Journal of Marine Science* 71:2141–2151.

Lin, Z.-S., and B.-L. Li. 2002. The maximum sustainable yield of Allee dynamic
system. *Ecological Modelling* 154:1–7.

McAllister, M. K., E. K. Pikitch, A. E. Punt, and R. Hilborn. 1994. A Bayesian
approach to stock assessment and harvest decisions using the
sampling/importance resampling algorithm. *Canadian Journal of Fisheries and
Aquatic Sciences* 51:2673–2687.

Romero, M. A., M. A. Coscarella, G. D. Adams, J. C. Pedraza, R. A. González, and
E. A. Crespo. 2022. Historical reconstruction of the population dynamics of
southern right whales in the southwestern Atlantic Ocean. *Scientific Reports*
12:3324.

Zerbini, A. N., E. J. Ward, P. G. Kinas, M. H. Engel, and A. Andriolo. 2011. A
Bayesian assessment of the conservation status of humpback whales (*Megaptera
novaeangliae*) in the western South Atlantic Ocean. *Journal of Cetacean Research
and Management* (Special Issue 3):131–144.

Zerbini, A. N., G. Adams, J. Best, P. J. Clapham, J. A. Jackson, and A. E. Punt.
2019. Assessing the recovery of an Antarctic predator from historical
exploitation. *Royal Society Open Science* 6(10):190368.
