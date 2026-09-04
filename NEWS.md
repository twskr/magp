# magp 0.12.0.9000

- Added a pkgdown website with a getting-started guide, an optimization
  workflow, a reproducible covariance-engine benchmark, and a curated function
  reference.
- Added package-level citation metadata for both the software and the
  associated methodology paper.
- Added `magp_bayes_optimize_from_scratch()` to connect simulated-annealing
  initial design, initial objective evaluation, and sequential MaGP Bayesian
  optimization in one call.
- Added reproducible handling of stochastic objective functions while
  preserving the caller's random-number state.
- Retained the complete initial design and its responses in the optimization
  result for review and reuse.
- Added early validation for design, fitting, acquisition, and objective
  controls before the initial-design search begins.

# magp 0.11.0.9000

- Added expected improvement for both two-dimensional and full-mapping models.
- Added mixed acquisition search across quantitative bounds and sequence
  permutations, with multi-start and optional parallel execution.
- Added a sequential Bayesian optimization interface for expensive objective
  functions.
- Added deterministic sequence sampling for cases where complete permutation
  enumeration would be too large.
- Added duplicate-input protection, early stopping, complete search history,
  and final-model refitting.

# magp 0.10.0.9000

- Added multi-start optimization to `magp2d_fit()` and `magpfull_fit()`.
- Added optional local parallel execution through portable socket clusters.
- Made seeded sequential and parallel fits use the same parameter starts.
- Added start-level objective values, convergence results, warnings, and errors
  to each fitted object.

# magp 0.9.0.9000

- Completed the quantitative-sequence initial-design workflow.
- Added `magp_quantitative_criterion()` and `magp_quantitative_design()` for
  maximin-style Latin hypercubes.
- Added `magp_joint_criterion()` to evaluate the combined separation of
  quantitative and sequence inputs.
- Added `magp_initial_design()` to construct and align both design portions
  while preserving their individual structures.
- Added reproducible simulated-annealing searches and validation across
  several run sizes and component counts.

# magp 0.8.0.9000

- Added tools for constructing the sequence portion of a quantitative-sequence
  initial design.
- Added `magp_sequence_criterion()` to assess pair balance and Hamming-distance
  space filling in a sequence design.
- Added `magp_sequence_design()` to optimize the sequence portion with
  reproducible simulated annealing.
- Added reproducible simulated-annealing search with configurable criterion
  weights and search settings.

# magp 0.8.0

- Prepared the two-dimensional and full-mapping models for the first CRAN
  release.
- Added a consistent `predict()` interface, including plug-in predictive
  uncertainty, for both fitted model classes.
- Moved the main covariance, gradient, and cross-covariance calculations to
  C++ with `Rcpp`.
- Added automatic input checks and reusable scaling for quantitative inputs.
- Set the default nugget variance to `tau = 0.001`.

# magp 0.7.0

- Renamed the package to `magp`.
- Added two-dimensional and full-mapping MaGP models.
- Added matching `predict()` methods for both fitted model classes.
- Used `Rcpp` for covariance, gradient, and cross-covariance calculations.
- Added automatic input-role detection, sequence-permutation checks, and
  stored min-max scaling for quantitative inputs outside `[0, 1]`.
