# Continue Bayesian optimization from completed experiments

Use this function when initial experiments and their responses are
already available. It fits a MaGP model, selects a new input with
expected improvement, evaluates `FUN`, and adds the new result to the
data. This process repeats for at most `n_iter` new evaluations.

## Usage

``` r
magp_bayes_optimize(
  FUN,
  X,
  y = NULL,
  model = c("2d", "full"),
  direction = c("minimize", "maximize"),
  n_iter = 10L,
  xi = 0,
  stop_ei = 0,
  stop_patience = 3L,
  seed = NULL,
  fit_control = list(),
  acquisition_control = list(),
  objective_args = list(),
  verbose = TRUE
)
```

## Arguments

- FUN:

  Function that evaluates one experiment. Its argument names must match
  the columns of `X`. It must return one finite number, or a list with
  one finite number named `Score` or `Value`.

- X:

  Initial inputs as a numeric matrix or data frame. The first half of
  the columns contains quantitative values. The second half contains the
  sequence positions, with each row forming a permutation of `1:q`. `X`
  may also contain a response column named `y`.

- y:

  Numeric responses for the rows of `X`. Leave this as `NULL` when `X`
  contains a column named `y`.

- model:

  MaGP mapping to fit. Use `"2d"` for the compact mapping or `"full"`
  for the full mapping.

- direction:

  Use `"minimize"` when smaller responses are better and `"maximize"`
  when larger responses are better.

- n_iter:

  Maximum number of new experiments to evaluate.

- xi:

  Nonnegative expected-improvement offset. The default, `0`, uses the
  current best response as the improvement target. Larger values require
  a candidate to exceed that target by more.

- stop_ei:

  Nonnegative early-stopping threshold for expected improvement.

- stop_patience:

  Number of consecutive selected points with expected improvement less
  than or equal to `stop_ei` required to stop early.

- seed:

  Optional nonnegative whole-number seed for reproducible model starts
  and acquisition searches. The caller's random-number state is
  preserved.

- fit_control:

  Optional named list passed to
  [`magp2d_fit()`](https://twskr.github.io/magp/reference/magp2d_fit.md)
  or
  [`magpfull_fit()`](https://twskr.github.io/magp/reference/magpfull_fit.md).
  Common choices include `tau`, `maxeval`, `n_starts`, and `workers`.
  This function supplies `X`, `y`, and `seed`.

- acquisition_control:

  Optional named list passed to
  [`magp_next_point()`](https://twskr.github.io/magp/reference/magp_next_point.md).
  Common choices include `lower`, `upper`, `sequences`, `n_starts`,
  `workers`, and `maxit`. This function supplies the fitted model,
  direction, `xi`, current best response, and seed.

- objective_args:

  Optional named list of fixed arguments passed to `FUN` in addition to
  the input columns.

- verbose:

  If `TRUE`, print the observed response and expected improvement after
  each new evaluation.

## Value

An object of class `magp_bayes_opt`. The final fitted model includes
every completed evaluation.

## How the loop works

The initial rows of `X` are treated as completed experiments and are not
evaluated again. Each iteration performs four steps:

1.  fit the selected MaGP model to all results collected so far;

2.  use
    [`magp_next_point()`](https://twskr.github.io/magp/reference/magp_next_point.md)
    to select an unobserved input;

3.  call `FUN` at that input; and

4.  add the response and refit the model.

`FUN` receives one named argument for each input column. For columns
`A`, `B`, `a`, and `b`, for example, the call is equivalent to
`FUN(A = ..., B = ..., a = ..., b = ...)`.

## Reading the result

The returned object contains:

- `call`: the function call;

- `best_point`: the input row with the best observed response;

- `best_value`: the best observed response;

- `best_index`: the row number of the best result in `X` and `y`;

- `history`: the initial and newly evaluated rows in evaluation order;

- `model`: the final fitted MaGP model;

- `X` and `y`: all inputs and responses used by the final model;

- `acquisitions`: details from each call to
  [`magp_next_point()`](https://twskr.github.io/magp/reference/magp_next_point.md);

- `mapping` and `direction`: the model and optimization direction;

- `iterations_requested` and `iterations_completed`: the requested and
  completed numbers of new evaluations;

- `stop_reason`: why the loop ended;

- `xi`, `stop_ei`, and `stop_patience`: the acquisition and stopping
  settings; and

- `fit_control` and `acquisition_control`: the control lists used in the
  search.

## References

Jones, D. R., Schonlau, M., and Welch, W. J. (1998). Efficient Global
Optimization of Expensive Black-Box Functions. Journal of Global
Optimization, 13, 455-492.
[doi:10.1023/A:1008306431147](https://doi.org/10.1023/A%3A1008306431147)
.

## Examples

``` r
# \donttest{
design <- magp_initial_design(n = 6, q = 3, seed = 1)
objective <- function(quantity_1, quantity_2, quantity_3,
                      sequence_1, sequence_2, sequence_3) {
  quantities <- c(quantity_1, quantity_2, quantity_3)
  positions <- c(sequence_1, sequence_2, sequence_3)
  -sum((quantities - c(0.2, 0.6, 0.8))^2) -
    0.02 * sum((positions - c(1, 3, 2))^2)
}
initial_y <- apply(design$design, 1L, function(row) {
  do.call(objective, as.list(row))
})
result <- magp_bayes_optimize(
  FUN = objective,
  X = design$design,
  y = initial_y,
  direction = "maximize",
  n_iter = 1,
  seed = 2,
  fit_control = list(maxeval = 100),
  acquisition_control = list(n_starts = 2, maxit = 20),
  verbose = FALSE
)
result$best_point
#> quantity_1 quantity_2 quantity_3 sequence_1 sequence_2 sequence_3 
#> 0.08333333 0.58333333 0.91666667 1.00000000 2.00000000 3.00000000 
result$best_value
#> [1] -0.0675
result$history
#>   Iteration Initial quantity_1 quantity_2 quantity_3 sequence_1 sequence_2
#> 1         0    TRUE 0.08333333 0.58333333 0.91666667          1          2
#> 2         0    TRUE 0.75000000 0.08333333 0.75000000          3          1
#> 3         0    TRUE 0.41666667 0.75000000 0.08333333          2          1
#> 4         0    TRUE 0.58333333 0.91666667 0.58333333          1          3
#> 5         0    TRUE 0.25000000 0.25000000 0.41666667          2          3
#> 6         0    TRUE 0.91666667 0.41666667 0.25000000          3          2
#> 7         1   FALSE 0.00000000 1.00000000 1.00000000          2          3
#>   sequence_3      Value ExpectedImprovement PredictedMean
#> 1          3 -0.0675000                  NA            NA
#> 2          2 -0.7319444                  NA            NA
#> 3          3 -0.7030556                  NA            NA
#> 4          2 -0.2941667                  NA            NA
#> 5          1 -0.3119444                  NA            NA
#> 6          1 -0.9697222                  NA            NA
#> 7          1 -0.2800000            1.793363    -0.4690127
#>   PredictedStandardError
#> 1                     NA
#> 2                     NA
#> 3                     NA
#> 4                     NA
#> 5                     NA
#> 6                     NA
#> 7               4.982345
# }
```
