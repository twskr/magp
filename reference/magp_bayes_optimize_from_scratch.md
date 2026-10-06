# Start Bayesian optimization before any experiments have been run

Use this function when no initial data are available. It creates an
initial quantitative-sequence design, evaluates `FUN` at those runs, and
then uses MaGP expected improvement to select later runs.

## Usage

``` r
magp_bayes_optimize_from_scratch(
  FUN,
  n_initial,
  q,
  model = c("2d", "full"),
  direction = c("minimize", "maximize"),
  n_iter = 10L,
  xi = 0,
  stop_ei = 0,
  stop_patience = 3L,
  seed = NULL,
  design_control = list(),
  fit_control = list(),
  acquisition_control = list(),
  objective_args = list(),
  verbose = TRUE
)
```

## Arguments

- FUN:

  Function that evaluates one experiment. It must accept the generated
  arguments `quantity_1`, ..., `quantity_q`, `sequence_1`, ...,
  `sequence_q`. It must return one finite number, or a list with one
  finite number named `Score` or `Value`.

- n_initial:

  Number of experiments in the generated initial design.

- q:

  Number of components. It must be at least three.

- model:

  MaGP mapping to fit. Use `"2d"` for the compact mapping or `"full"`
  for the full mapping.

- direction:

  Use `"minimize"` when smaller responses are better and `"maximize"`
  when larger responses are better.

- n_iter:

  Maximum number of new experiments selected after the initial design
  has been evaluated.

- xi:

  Nonnegative expected-improvement offset. The default, `0`, uses the
  current best response as the improvement target.

- stop_ei:

  Nonnegative early-stopping threshold for expected improvement.

- stop_patience:

  Number of consecutive selected points with expected improvement less
  than or equal to `stop_ei` required to stop early.

- seed:

  Optional nonnegative whole-number seed. It makes the initial design,
  model starts, acquisition searches, and random values produced by
  `FUN` reproducible while preserving the caller's random-number state.

- design_control:

  Optional named list passed to
  [`magp_initial_design()`](https://twskr.github.io/magp/reference/magp_initial_design.md).
  Use `sequence_method` to choose `"random"`, `"sfta"`, or `"sann"`.
  Other common choices include `sequence_maxit`, `sfta_control`,
  `quantity_maxit`, and `alignment_maxit`. This function supplies `n`,
  `q`, and `seed`.

- fit_control:

  Optional named list passed to
  [`magp2d_fit()`](https://twskr.github.io/magp/reference/magp2d_fit.md)
  or
  [`magpfull_fit()`](https://twskr.github.io/magp/reference/magpfull_fit.md).
  Common choices are `maxeval`, `n_starts`, and `workers`. This function
  supplies `X`, `y`, and `seed`.

- acquisition_control:

  Optional named list passed to
  [`magp_next_point()`](https://twskr.github.io/magp/reference/magp_next_point.md).
  Common choices are `lower`, `upper`, `sequences`, `n_starts`,
  `workers`, and `maxit`. This function supplies the fitted model,
  direction, `xi`, current best response, and seed.

- objective_args:

  Optional named list of fixed arguments passed to `FUN` in addition to
  the generated input columns.

- verbose:

  If `TRUE`, print one line after each initial and sequential
  evaluation.

## Value

An object of class `magp_bayes_opt` containing the initial design, all
completed evaluations, and the final fitted model.

## How to use this function

Define `FUN`, choose `n_initial` and `q`, and state whether the response
should be minimized or maximized. The function then:

1.  creates an initial design;

2.  calls `FUN` once for each initial row;

3.  fits the selected MaGP model; and

4.  selects and evaluates as many as `n_iter` additional rows.

The generated columns are named `quantity_1` through `quantity_q` and
`sequence_1` through `sequence_q`. Use
[`magp_bayes_optimize()`](https://twskr.github.io/magp/reference/magp_bayes_optimize.md)
instead when initial experiments and responses already exist.

## Reading the result

The result has the same main components as
[`magp_bayes_optimize()`](https://twskr.github.io/magp/reference/magp_bayes_optimize.md),
including `best_point`, `best_value`, `history`, the final `model`, and
all recorded search settings. It also contains:

- `initial_design`: the generated quantitative-sequence design;

- `initial_response`: the response from each initial run;

- `initial_evaluations`: the number of initial runs; and

- `design_control`: the initial-design settings that were used; and

- `started_from_initial_design`: `TRUE`, marking that the workflow
  generated its own starting design.

## Examples

``` r
# \donttest{
objective <- function(quantity_1, quantity_2, quantity_3,
                      sequence_1, sequence_2, sequence_3) {
  quantities <- c(quantity_1, quantity_2, quantity_3)
  sequence <- c(sequence_1, sequence_2, sequence_3)
  -sum((quantities - c(0.2, 0.6, 0.8))^2) -
    0.01 * sum((sequence - c(1, 3, 2))^2)
}
result <- magp_bayes_optimize_from_scratch(
  FUN = objective,
  n_initial = 6,
  q = 3,
  direction = "maximize",
  n_iter = 1,
  seed = 1,
  design_control = list(
    sequence_maxit = 100,
    quantity_maxit = 100,
    alignment_maxit = 100
  ),
  fit_control = list(maxeval = 100),
  acquisition_control = list(n_starts = 2, maxit = 20),
  verbose = FALSE
)
result$initial_design$design
#>      quantity_1 quantity_2 quantity_3 sequence_1 sequence_2 sequence_3
#> [1,] 0.75000000 0.08333333 0.25000000          1          2          3
#> [2,] 0.08333333 0.25000000 0.41666667          3          1          2
#> [3,] 0.25000000 0.58333333 0.91666667          2          1          3
#> [4,] 0.91666667 0.41666667 0.75000000          1          3          2
#> [5,] 0.41666667 0.75000000 0.08333333          2          3          1
#> [6,] 0.58333333 0.91666667 0.58333333          3          2          1
result$best_point
#> quantity_1 quantity_2 quantity_3 sequence_1 sequence_2 sequence_3 
#>  0.2500000  0.5833333  0.9166667  2.0000000  1.0000000  3.0000000 
result$best_value
#> [1] -0.07638889
result$history
#>   Iteration Initial quantity_1 quantity_2 quantity_3 sequence_1 sequence_2
#> 1         0    TRUE 0.75000000 0.08333333 0.25000000          1          2
#> 2         0    TRUE 0.08333333 0.25000000 0.41666667          3          1
#> 3         0    TRUE 0.25000000 0.58333333 0.91666667          2          1
#> 4         0    TRUE 0.91666667 0.41666667 0.75000000          1          3
#> 5         0    TRUE 0.41666667 0.75000000 0.08333333          2          3
#> 6         0    TRUE 0.58333333 0.91666667 0.58333333          3          2
#> 7         1   FALSE 0.00000000 0.00000000 1.00000000          1          3
#>   sequence_3       Value ExpectedImprovement PredictedMean
#> 1          3 -0.89194444                  NA            NA
#> 2          2 -0.36305556                  NA            NA
#> 3          3 -0.07638889                  NA            NA
#> 4          2 -0.54972222                  NA            NA
#> 5          1 -0.60305556                  NA            NA
#> 6          1 -0.35416667                  NA            NA
#> 7          2 -0.44000000            2.165427    -0.4700069
#>   PredictedStandardError
#> 1                     NA
#> 2                     NA
#> 3                     NA
#> 4                     NA
#> 5                     NA
#> 6                     NA
#> 7               5.908141
# }
```
