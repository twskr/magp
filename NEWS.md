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
