# Run sequential Bayesian optimization with a MaGP surrogate

Repeatedly fits a two-dimensional or full-mapping MaGP model, maximizes
expected improvement over quantitative inputs and sequence permutations,
evaluates a user-supplied objective function, and adds the result to the
training data.

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

  Objective function. It is called with one named argument for each
  input column. It must return one finite numeric value or a list
  containing a finite numeric `Score` or `Value` component.

- X:

  Initial quantitative-sequence inputs. A response column named `y` may
  be included.

- y:

  Optional initial response vector. It may be omitted when `X` contains
  a column named `y`.

- model:

  Sequence mapping used by the surrogate: `"2d"` or `"full"`.

- direction:

  Whether `FUN` is being minimized or maximized.

- n_iter:

  Maximum number of new objective evaluations.

- xi:

  Nonnegative exploration offset for expected improvement.

- stop_ei:

  Nonnegative absolute expected-improvement threshold. The search stops
  after `stop_patience` consecutive evaluated points at or below this
  threshold.

- stop_patience:

  Positive whole number of consecutive low-EI evaluations required for
  early stopping.

- seed:

  Optional nonnegative whole-number seed. It controls model starts,
  sequence sampling, and acquisition starts without changing the
  caller's random-number state.

- fit_control:

  Named list of additional arguments for
  [`magp2d_fit()`](https://twskr.github.io/magp/reference/magp2d_fit.md)
  or
  [`magpfull_fit()`](https://twskr.github.io/magp/reference/magpfull_fit.md),
  such as `tau`, `n_starts`, and `workers`. `X`, `y`, and `seed` are
  managed by this function.

- acquisition_control:

  Named list of additional arguments for
  [`magp_next_point()`](https://twskr.github.io/magp/reference/magp_next_point.md),
  such as `lower`, `upper`, `sequences`, `n_starts`, `workers`, and
  `maxit`. The model, direction, `xi`, reference value, and seed are
  managed by this function.

- objective_args:

  Named list of fixed additional arguments supplied to `FUN` after the
  input columns.

- verbose:

  Logical; if `TRUE`, print one line after every new objective
  evaluation.

## Value

An object of class `magp_bayes_opt` with components `best_point`,
`best_value`, `history`, `model`, `X`, `y`, `acquisitions`, and the
search settings. The final fitted model includes every completed
evaluation.

## Details

The initial responses are treated as completed experiments and are not
evaluated again. At each iteration, the current observed optimum is used
as the reference value for expected improvement. Previously observed
inputs are excluded by default through
[`magp_next_point()`](https://twskr.github.io/magp/reference/magp_next_point.md).

`FUN` follows a named-argument interface. For input columns `A`, `B`,
`a`, and `b`, for instance, the function is called as
`FUN(A, B, a, b, ...)`. A database lookup or other data source can be
used by wrapping it in a function with the same interface.

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
  objective,
  design$design,
  initial_y,
  direction = "maximize",
  n_iter = 1,
  seed = 2,
  fit_control = list(maxeval = 50),
  acquisition_control = list(n_starts = 2, maxit = 20)
)
#> Warning: nloptr did not report convergence; status 5: NLOPT_MAXEVAL_REACHED: Optimization stopped because maxeval (above) was reached.. Inspect the fitted model before using it.
#> Warning: nloptr did not report convergence; status 5: NLOPT_MAXEVAL_REACHED: Optimization stopped because maxeval (above) was reached.. Inspect the fitted model before using it.
#> iteration 1 value -1 EI 2.08444 
result$best_point
#> quantity_1 quantity_2 quantity_3 sequence_1 sequence_2 sequence_3 
#> 0.08333333 0.58333333 0.91666667 1.00000000 2.00000000 3.00000000 
result$best_value
#> [1] -0.0675
# }
```
