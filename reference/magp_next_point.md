# Select the next quantitative-sequence experiment

Searches the allowed quantitative values and sequence permutations, then
returns the unobserved point with the largest expected improvement.

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

  Use `"minimize"` when smaller responses are better and `"maximize"`
  when larger responses are better.

- xi:

  Nonnegative improvement offset used in expected improvement.

- best:

  Optional response that a new point should improve upon. When it is
  omitted, the function uses the best observed response in `object`.

- lower, upper:

  Optional lower and upper bounds for the quantitative inputs. Supply
  one value for all components or one value per component. The defaults
  use the prediction ranges stored in the fitted model.

- sequences:

  Optional matrix of candidate sequence permutations. When omitted,
  every permutation is used if their number does not exceed
  `max_sequences`; otherwise a reproducible sample is searched.

- max_sequences:

  Largest number of sequence candidates generated when `sequences` is
  not supplied.

- n_starts:

  Number of quantitative starting points searched for each sequence
  candidate.

- workers:

  Number of local worker processes. Use `1` for sequential execution. At
  most two processes are used.

- maxit:

  Maximum optimization iterations for each quantitative start.

- factr, pgtol:

  Advanced convergence settings passed to
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html) for its
  `"L-BFGS-B"` method.

- exclude_observed:

  If `TRUE`, do not return an input that is already in the training
  data.

- duplicate_tolerance:

  Nonnegative absolute tolerance used to identify an observed input
  after the model's quantitative scaling is applied.

- seed:

  Optional nonnegative whole-number seed for sequence sampling and
  quantitative starting points.

## Value

An object of class `magp_next_point` containing the selected point, its
prediction, and search diagnostics.

## Sequence search

If `sequences` is supplied, only those rows are searched. Otherwise, the
function searches every permutation when there are no more than
`max_sequences`. For a larger sequence space, it searches a reproducible
sample of `max_sequences` permutations when `seed` is supplied.

## Reading the result

The returned object contains:

- `call`: the function call;

- `point`: the selected input as a one-row data frame;

- `expected_improvement`: the score of the selected point;

- `predicted_mean` and `predicted_standard_error`: the MaGP prediction;

- `direction`, `best_observed`, and `xi`: the expected-improvement
  settings;

- `bounds`: the quantitative lower and upper bounds;

- `sequences` and `sequence_source`: the permutations searched and how
  they were obtained;

- `selected_sequence` and `selected_start`: the winning search indices;

- `n_starts`: the number of quantitative starts per sequence;

- `workers_requested`, `workers_used`, and `execution`: the
  parallel-search settings actually used; and

- `diagnostics`: the result of every sequence and starting-point search.

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
next_run$expected_improvement
#> [1] 7.521459
# }
```
