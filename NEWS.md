# SAEforestArrow 2.0.0

Scalable fork of SAEforest 1.0.1 for very large (raster) populations.

## New

* `SAEpopulation()`: population from a raster (all layers or a subset as covariates, a mask
  layer, domains from polygons rasterized by cell centre). Valid cells are extracted once,
  in parallel and block-wise, into an on-disk `fst` cache that is streamed in chunks.
  Also wraps in-memory data.frames and re-opens existing caches.
* `SAEforest_model()` accepts `SAEpop` objects or rasters (`pop_options`) as `pop_data`,
  and gains `total_factor` (domain totals and their MSE), `cores`, `chunk_rows`,
  `forest_batch`, `seed`, `verbose` and `engine`.
* Scalable engine for domain means (`meanOnly = TRUE`, unit-level data), used by default:
  memory no longer grows with population size x `B`. Point estimates are identical to the
  original implementation; the nonparametric MSE bootstrap follows the same distribution.
  - only domain sums of forest predictions are kept, the population is streamed;
  - random effects via `ranef()` look-up (identical to `predict.merMod` on new data);
  - bootstrap unit errors are drawn in C++ one replicate at a time and summed per domain;
  - bootstrap MERFs are fitted in parallel on the (small) bootstrap samples only.
* `flatten_forest()`, `predict_flat_forest()`: multi-threaded C++ prediction of ranger
  regression forests, bit-identical to `ranger::predict()` and ~10-15x faster with many
  covariates (ranger sorts all covariates at each prediction call).
* `summarize_indicators()` reports `Total` besides `Mean` when available.
* `SAEaggregate()`: means and totals of any aggregation of domains (e.g. municipalities to
  provinces, regions, country) with exact bootstrap MSEs, computed from the bootstrap
  replicates stored in the model (covariances between domains included). It is identical to
  re-running the model with `aggregate_to`.
* Box-Cox transformed MERF (scalable engine): arguments `transform = "boxcox"`, `bc_lambda`,
  `bc_shift`, `bt_nodes`, `bt_grid`.
  - The MERF is fitted to the transformed target.
  - Domain means are empirical-best-type predictions on the original scale: expected
    back-transformed values given the predicted random effect, integrating the unit error
    and the remaining random-effect uncertainty by Gauss-Hermite quadrature (Molina & Rao,
    2010).
  - The nonparametric bootstrap runs on the transformed scale; true and estimated
    bootstrap domain means are back-transformed unit by unit.
  - With `bc_lambda = 1`, `bc_shift = 1` (identity) the results equal `transform = "none"`
    with the same bootstrap error pool, as long as no transformed value falls below
    `-bc_shift`, where the back-transformation is truncated.
* `boot_error_pool` (`"auto"`, `"domain"`, `"all"`): pool of the unit errors of the
  bootstrap populations.
  - `"domain"` (the original scheme) draws the errors of sampled domains from their own,
    within-domain centred residuals; `"all"` draws them from all residuals.
  - `"auto"` uses `"domain"` without a transformation (results unchanged) and `"all"` with
    Box-Cox. After a back-transformation the true domain means depend on the whole error
    distribution, which the few own residuals of a sampled domain do not represent.
* `boot_hetero` (`"none"`, `"unit"`, `"unit_re"`): optional heteroscedastic variants of the
  bootstrap (unit errors with SD proportional to a power of the level; optionally random
  effects proportional to the domain level). `"none"` (default) is the original bootstrap.

## Differences from SAEforest 1.0.1

* `MERFranger()`: a unit that is in-bag in every tree (possible with few trees) has no
  out-of-bag prediction. `ranger` returns `NaN`, `lmer` then drops the unit, and the fit stopped
  with "Missing data in dependent variable". Such units now get the prediction of the whole
  forest. Results are unchanged whenever all out-of-bag predictions exist (always in practice
  with a few hundred trees).

* Block bootstrap of errors: in the original code, a domain with a single sample unit
  calls `sample(x, ...)` with a length-one numeric `x`, so that its errors are drawn from
  `1:x` instead of `x`. The scalable engine draws the errors in C++ and is not affected.
  `engine = "legacy"` keeps the original code unchanged, including this behaviour, for
  comparison.

# SAEforest 1.0.1

# - remove dependencies to rgeos, maptools and sp
