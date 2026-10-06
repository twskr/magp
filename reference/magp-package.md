# magp: Mapping-Based Additive Gaussian Process Models

Fits additive Gaussian process models for experiments in which each
component has both a quantitative level and a position in a sequence.

## Details

For `q` components, the first `q` input columns contain quantitative
levels and the next `q` columns contain sequence positions. Every
sequence row must be a permutation of `1:q`. Quantitative columns
outside `[0, 1]` are scaled during fitting, and the same transformation
is used for new data.

[`magp2d_fit()`](https://twskr.github.io/magp/reference/magp2d_fit.md)
represents the sequence positions in two latent dimensions.
[`magpfull_fit()`](https://twskr.github.io/magp/reference/magpfull_fit.md)
uses `q - 1` latent dimensions. Both fitted model classes support
[`stats::predict()`](https://rdrr.io/r/stats/predict.html) for point
predictions and plug-in predictive uncertainty. The computational
kernels for covariance matrices, analytical gradients, and
cross-covariances are implemented in C++ with `Rcpp`.

[`magp_initial_design()`](https://twskr.github.io/magp/reference/magp_initial_design.md)
constructs a space-filling quantitative-sequence design before responses
are collected. It combines a Latin hypercube with sequence permutations
generated randomly or improved with simulated annealing or space-filling
threshold accepting. Joint alignment preserves both component designs.

[`magp_expected_improvement()`](https://twskr.github.io/magp/reference/magp_expected_improvement.md)
evaluates improvement using latent predictive uncertainty.
[`magp_next_point()`](https://twskr.github.io/magp/reference/magp_next_point.md)
searches quantitative bounds and sequence permutations for the next
experiment, while
[`magp_bayes_optimize()`](https://twskr.github.io/magp/reference/magp_bayes_optimize.md)
runs the sequential fitting and evaluation loop from completed
experiments.
[`magp_bayes_optimize_from_scratch()`](https://twskr.github.io/magp/reference/magp_bayes_optimize_from_scratch.md)
first generates and evaluates an initial design, then continues through
the same sequential loop.

## See also

Useful links:

- <https://github.com/twskr/magp>

- <https://CRAN.R-project.org/package=magp>

- Report bugs at <https://github.com/twskr/magp/issues>

## Author

**Maintainer**: Tony Wang <wangtony883@gmail.com> \[copyright holder\]

Authors:

- Tony Wang <wangtony883@gmail.com> \[copyright holder\]

- Qian Xiao \[copyright holder\]

Other contributors:

- Yaping Wang \[copyright holder\]

- Abhyuday Mandal \[copyright holder\]

- Xinwei Deng \[copyright holder\]
