# Find the next quantitative-sequence experiment

Maximizes expected improvement over quantitative inputs and sequence
permutations. Each sequence is searched from several quantitative
starting points, and the best result across all searches is returned.

## Usage

``` r
magp_next_point(
  object,
  direction = c("minimize", "maximize"),
  xi = 0,
  best = NULL,
  lower = NULL,
  upper = NULL,
  sequences = NULL,
  max_sequences = 120L,
  n_starts = 5L,
  workers = 1L,
  maxit = 100L,
  factr = 1e+07,
  pgtol = 0,
  exclude_observed = TRUE,
  duplicate_tolerance = sqrt(.Machine$double.eps),
  seed = NULL
)
```

## Arguments

- object:

  A fitted `magp2d` or `magpfull` model.

- direction:

  Whether the objective is being minimized or maximized.

- xi:

  A nonnegative exploration offset used by expected improvement.

- best:

  Optional finite reference value. The current observed optimum is used
  when omitted.

- lower, upper:

  Optional quantitative bounds. Each may contain one value or `q`
  values. Defaults are `[0, 1]` for inputs fitted on that scale and the
  stored training range for inputs that were scaled during fitting.

- sequences:

  Optional matrix of candidate sequence permutations. When omitted,
  every permutation is used if their number does not exceed
  `max_sequences`; otherwise a reproducible sample is searched.

- max_sequences:

  Maximum number of automatically generated sequence candidates.

- n_starts:

  Number of quantitative starting points used for each sequence. The
  first is the midpoint of the bounds and the rest are random.

- workers:

  Number of local worker processes. Values greater than one use a
  portable socket cluster and are capped at the number of search tasks.

- maxit:

  Maximum number of `L-BFGS-B` iterations for each start.

- factr, pgtol:

  Convergence controls passed to
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html) for its
  `L-BFGS-B` method.

- exclude_observed:

  Logical; if `TRUE`, previously observed inputs are not eligible for
  selection.

- duplicate_tolerance:

  Nonnegative absolute tolerance used to identify an observed input
  after the model's quantitative scaling is applied.

- seed:

  Optional nonnegative whole-number seed for sequence sampling and
  quantitative starting points.

## Value

An object of class `magp_next_point`. Its `point` component is a one-row
data frame ready for evaluation. The object also contains expected
improvement, predictive mean and standard error, the searched sequences,
and start-level diagnostics.

## Examples

``` r
# \donttest{
train <- read.table(
  system.file("extdata", "example_train.txt", package = "magp"),
  header = TRUE
)
fit <- magp2d_fit(train, seed = 1)
next_run <- magp_next_point(
  fit,
  direction = "minimize",
  n_starts = 2,
  maxit = 20,
  seed = 2
)
next_run$point
#>          A       B         C         D a b c d
#> 1 0.678576 0.98949 0.5597923 0.5015999 1 4 3 2
# }
```
