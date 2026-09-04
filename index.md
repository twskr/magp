# magp

[![CRAN
status](https://www.r-pkg.org/badges/version/magp)](https://CRAN.R-project.org/package=magp)
[![R-CMD-check](https://github.com/twskr/magp/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/twskr/magp/actions/workflows/R-CMD-check.yaml)
[![Documentation](https://img.shields.io/badge/docs-pkgdown-2c7fb8.svg)](https://twskr.github.io/magp/)

The stable release is available from CRAN. This branch contains the
development version, `0.12.0.9000`.

`magp` is designed for experiments in which every component has both an
amount and a position in a sequence. It fits an additive Gaussian
process that uses both parts of the input, then returns predictions with
optional uncertainty estimates. You can choose a compact two-dimensional
mapping or a full mapping with `q - 1` dimensions. The most intensive
covariance and gradient calculations run in C++ through `Rcpp`.

## Installation

Install the stable release from CRAN:

``` r

install.packages("magp")
library(magp)
```

To install the current development branch from GitHub:

``` r

install.packages("remotes")
remotes::install_github("twskr/magp")
```

A C++ toolchain is required when installing from source or GitHub.

The complete function reference and tutorials are available on the
[package website](https://twskr.github.io/magp/).

## Citation

Run `citation("magp")` for the software citation and the methodology
paper. The stable package has the CRAN DOI
[10.32614/CRAN.package.magp](https://doi.org/10.32614/CRAN.package.magp).

## Data format

For `q` components, the input must contain `2*q` columns:

1.  the first `q` columns contain quantitative inputs;
2.  the last `q` columns contain sequence positions.

Every row in the sequence columns must contain each value from `1` to
`q` exactly once. A response column named `y` may be included in the
same data frame; when it is present, the fitting functions can identify
both `y` and `q` automatically.

Quantitative columns outside `[0, 1]` are transformed by min-max scaling
during fitting. The fitted ranges are saved and used again for
prediction. Inputs already in `[0, 1]` are not changed.

## Example

``` r

train <- read.table(
  system.file("extdata", "example_train.txt", package = "magp"),
  header = TRUE
)
test <- read.table(
  system.file("extdata", "example_test.txt", package = "magp"),
  header = TRUE
)

fit_2d <- magp2d_fit(train, tau = 0.001, seed = 1)
prediction_2d <- predict(fit_2d, test)
magp2d_rmse(prediction_2d, test$y)

fit_full <- magpfull_fit(train, tau = 0.001, seed = 1)
prediction_full <- predict(fit_full, test)
magp2d_rmse(prediction_full, test$y)

uncertainty <- predict(
  fit_2d,
  test[1:5, ],
  se.fit = TRUE,
  type = "response"
)
data.frame(
  prediction = uncertainty$fit,
  standard_error = uncertainty$se.fit
)
```

`tau` is a fixed nugget variance added to the covariance diagonal. It is
a variance, not a standard deviation.

## Prediction types

- `"script"` reproduces the fitted training-row convention when a new
  row exactly matches a training row.
- `"response"` includes the nugget when calculating uncertainty for a
  future response.
- `"latent"` returns uncertainty for the noise-free surface.

The reported standard errors treat the fitted covariance parameters as
fixed.

## Multi-start fitting

Both models can be fitted from several initial parameter vectors. This
is useful when one optimizer run may settle at a local solution. Set
`workers` above one to run the starts in separate local R processes.

``` r

fit <- magp2d_fit(
  train,
  seed = 1,
  n_starts = 8,
  workers = 4
)
fit$multistart$starts
```

The returned model is the converged start with the lowest objective
value. If none of the starts converges, the lowest finite result is
returned with a warning. The start table records every objective value,
convergence code, warning, and error. A supplied seed gives the same set
of starts in sequential and parallel runs. `workers = 1` uses ordinary
sequential execution; larger values use a local socket cluster and are
capped at `n_starts`.

## Initial designs

The package can construct a complete quantitative-sequence initial
design. It first improves the sequence permutations, then generates a
space-filling Latin hypercube, and finally pairs the rows of the two
portions. The final step changes only the row order of the Latin
hypercube, so the individual quantity and sequence designs remain valid.

``` r

initial_design <- magp_initial_design(
  n = 16,
  q = 4,
  seed = 1
)
initial_design$design
initial_design$criteria
```

The first four columns in this example contain quantitative levels and
the last four contain sequence positions. The returned matrix can be
passed directly to either fitting function after a response vector has
been obtained.

The component criteria can also be used separately:

``` r

magp_sequence_criterion(initial_design$sequence)
magp_quantitative_criterion(initial_design$quantity)
magp_joint_criterion(
  initial_design$quantity,
  initial_design$sequence
)
```

All design criteria are minimized. A supplied seed makes each search
reproducible without changing the caller’s random-number state.

## Expected improvement and the next experiment

Expected improvement uses the latent predictive distribution from either
mapping model. Set `direction` explicitly so the observed optimum is
handled correctly.

``` r

ei <- magp_expected_improvement(
  fit_2d,
  test[1:5, ],
  direction = "minimize"
)

next_run <- magp_next_point(
  fit_2d,
  direction = "minimize",
  n_starts = 5,
  workers = 2,
  seed = 2
)
next_run$point
next_run$expected_improvement
```

For four components, the default search checks all 24 sequence
permutations. When the complete set is larger than `max_sequences`, the
function searches a reproducible sample instead. A specific set of
permutations can be supplied through `sequences`. Quantitative inputs
are optimized within the model’s stored prediction ranges, and completed
experiments are excluded by default.

## Sequential Bayesian optimization

[`magp_bayes_optimize()`](https://twskr.github.io/magp/reference/magp_bayes_optimize.md)
fits a MaGP surrogate, finds the point with the largest expected
improvement, evaluates the objective, and refits the model. The
objective receives one named argument for every input column and may
return a numeric value or a list containing `Score` or `Value`.

``` r

objective <- function(A, B, C, D, a, b, c, d) {
  quantities <- c(A, B, C, D)
  sequence <- c(a, b, c, d)
  target_quantities <- c(0.2, 0.4, 0.7, 0.9)
  target_sequence <- c(1, 3, 4, 2)
  -sum((quantities - target_quantities)^2) -
    0.01 * sum((sequence - target_sequence)^2)
}
initial_x <- train[, c("A", "B", "C", "D", "a", "b", "c", "d")]
initial_y <- apply(initial_x, 1L, function(row) {
  do.call(objective, as.list(row))
})

result <- magp_bayes_optimize(
  objective,
  initial_x,
  initial_y,
  direction = "maximize",
  n_iter = 3,
  seed = 3,
  fit_control = list(n_starts = 4, workers = 2),
  acquisition_control = list(n_starts = 5, workers = 2)
)
result$best_point
result$best_value
result$history
```

`stop_ei` is an absolute expected-improvement threshold. The loop stops
after `stop_patience` consecutive evaluated points at or below that
value. The returned model always includes every evaluation recorded in
`history`.

## Starting from a generated design

When no experiments have been run yet,
[`magp_bayes_optimize_from_scratch()`](https://twskr.github.io/magp/reference/magp_bayes_optimize_from_scratch.md)
connects the complete workflow. It builds the quantitative-sequence
initial design, evaluates the objective at those runs, and then
continues with sequential expected-improvement search.

``` r

objective <- function(quantity_1, quantity_2, quantity_3,
                      sequence_1, sequence_2, sequence_3) {
  quantities <- c(quantity_1, quantity_2, quantity_3)
  sequence <- c(sequence_1, sequence_2, sequence_3)
  -sum((quantities - c(0.2, 0.6, 0.8))^2) -
    0.01 * sum((sequence - c(1, 3, 2))^2)
}

result <- magp_bayes_optimize_from_scratch(
  objective,
  n_initial = 8,
  q = 3,
  direction = "maximize",
  n_iter = 3,
  seed = 4,
  design_control = list(
    sequence_maxit = 500,
    quantity_maxit = 500,
    alignment_maxit = 500
  ),
  fit_control = list(n_starts = 4, workers = 2),
  acquisition_control = list(n_starts = 5, workers = 2)
)
result$initial_design$design
result$initial_response
result$best_point
result$history
```

Use `design_control` to change the simulated-annealing settings for the
initial design. Use
[`magp_bayes_optimize()`](https://twskr.github.io/magp/reference/magp_bayes_optimize.md)
when initial experiments and responses are already available.
