# Design and sequential optimization

The development version of `magp` can connect three stages of an
experiment:

1.  generate an initial quantitative-sequence design;
2.  evaluate the response at those initial runs;
3.  use expected improvement to select later runs one at a time.

The objective function in this article is deterministic and inexpensive.
In a real study, the same function interface may call a simulator, run a
laboratory procedure through an external system, or return a measured
response entered after an experiment is complete.

## Build an initial design

``` r

library(magp)

design <- magp_initial_design(
  n = 8,
  q = 3,
  sequence_maxit = 200,
  quantity_maxit = 200,
  alignment_maxit = 200,
  seed = 4
)

design$design
#>      quantity_1 quantity_2 quantity_3 sequence_1 sequence_2 sequence_3
#> [1,]     0.4375     0.6875     0.9375          2          1          3
#> [2,]     0.0625     0.3125     0.6875          1          3          2
#> [3,]     0.6875     0.9375     0.4375          1          3          2
#> [4,]     0.3125     0.1875     0.0625          3          2          1
#> [5,]     0.9375     0.5625     0.8125          3          2          1
#> [6,]     0.5625     0.0625     0.5625          2          3          1
#> [7,]     0.1875     0.8125     0.3125          3          1          2
#> [8,]     0.8125     0.4375     0.1875          1          2          3
design$criteria
#>               sequence           quantitative joint_before_alignment 
#>              1.0318297              2.0769333              0.7997111 
#>                  joint 
#>              0.7052987
```

The first three columns contain the quantitative levels. The last three
columns contain sequence positions. Each quantitative column remains a
Latin hypercube, and every sequence row remains a valid permutation.

## Define an objective

The integrated interface passes one named argument for every input
column. This example rewards quantitative levels near a target and a
preferred sequence.

``` r

objective <- function(quantity_1, quantity_2, quantity_3,
                      sequence_1, sequence_2, sequence_3) {
  quantity <- c(quantity_1, quantity_2, quantity_3)
  sequence <- c(sequence_1, sequence_2, sequence_3)

  -sum((quantity - c(0.2, 0.6, 0.8))^2) -
    0.01 * sum((sequence - c(1, 3, 2))^2)
}
```

## Run the complete workflow

``` r

result <- magp_bayes_optimize_from_scratch(
  objective,
  n_initial = 8,
  q = 3,
  model = "2d",
  direction = "maximize",
  n_iter = 5,
  seed = 4,
  design_control = list(
    sequence_maxit = 500,
    quantity_maxit = 500,
    alignment_maxit = 500
  ),
  fit_control = list(
    n_starts = 4,
    workers = 2
  ),
  acquisition_control = list(
    n_starts = 5,
    workers = 2
  )
)

result$best_point
result$best_value
result$history
```

`direction` must state whether smaller or larger responses are
preferred. The returned history keeps the initial and sequential runs
together, so each selection can be reviewed and reproduced.

## Continue from completed experiments

When initial runs and responses already exist, use
[`magp_bayes_optimize()`](https://twskr.github.io/magp/reference/magp_bayes_optimize.md)
instead. It fits the selected MaGP surrogate, maximizes expected
improvement, evaluates the objective, and refits after each new result.

``` r

X_initial <- design$design
y_initial <- apply(X_initial, 1L, function(row) {
  do.call(objective, as.list(row))
})

result <- magp_bayes_optimize(
  objective,
  X = X_initial,
  y = y_initial,
  direction = "maximize",
  n_iter = 5,
  seed = 4
)
```

For laboratory work, do not connect the function directly to equipment
until the proposed point, bounds, units, and safety constraints have
been reviewed. The package selects statistical candidates; it does not
replace experimental oversight.
