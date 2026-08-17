# StateSpaceSIR 1.1.0

Version used for Romero et al. (in review, *Ecology*), *Accounting for
depensation in the population dynamics of southern right whales in the
southwestern Atlantic Ocean*.

## Depensation

* `StateSpaceSIR()` gains `allee_model`, selecting the depensation function
  applied to surplus production: `0` none, `1` Hilborn et al. (2014),
  `2` logistic, `3` Lin & Li (2002), `4` Haider et al. (2017).
* New parameter `P50` (`Pd` in the manuscript), the depletion level governing the
  strength of the Allee effect. Set a prior with `make_prior()` as for any other
  parameter.
* `Pmsy` is converted to the shape parameter `z` numerically, with a separate
  solver per depensation function (`pmsy_z_hilborn()`, `pmsy_z_logistic()`,
  `pmsy_z_linli()`, `pmsy_z_haider()`).
* `plot_production_function()` draws the realised surplus production curve.

## Model comparison

* `waic()` computes WAIC across a list of fitted models.
* Log-likelihoods are now reported per observation, so WAIC can be computed for
  the multivariate lognormal index likelihood.
* `bayes_factor()` and `weight_model()` handle the depensation models.

## Fixes

* Corrected the beta label in `summary_table()` and `plot_density()`. `"$\beta_"`
  was written without escaping the backslash, so `\b` was interpreted as a
  backspace and the catchability-trend parameter appeared as `$<BS>eta_{q_{flt1}}$`
  in summary tables and density-plot axes. Now `"$\\beta_"`.
* Fixed plotting of models with more than one index of abundance.
* Fixed a bounding issue in `run_SIR()` that could return abundances below the
  minimum-population constraint.
* Error rather than silent failure when a catchability drawn from the prior is
  `<= 0`.

# StateSpaceSIR 1.0.0

Version used for Romero et al. (2022),
[doi:10.1038/s41598-022-07370-6](https://doi.org/10.1038/s41598-022-07370-6).

* State-space dynamics with lognormal process error and an inverse-gamma prior
  on the process-error variance.
* Multiple low–high catch streams with a single estimated catch parameter, and
  period-specific struck-and-loss multipliers.
* Absolute abundance, relative abundance, and count-data likelihoods; catchability
  integrated out analytically, including the multivariate lognormal case.
* Prior on either `Pmsy` or `z`.
* Adapted from the HumpbackSIR package of Zerbini et al. (2011, 2019).
