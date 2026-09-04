# Construct the quantitative portion of an initial design

Uses simulated annealing to improve a Latin hypercube under a
maximin-style Euclidean-distance criterion. A candidate is created by
swapping two values within one column, which preserves the Latin
hypercube property.

## Usage

``` r
magp_quantitative_design(
  n,
  q,
  p = 15L,
  maxit = 10000L,
  temp = 0.01,
  tmax = 10L,
  initial = NULL,
  seed = NULL
)
```

## Arguments

- n:

  Number of design runs. Must be at least two.

- q:

  Number of quantitative variables. Must be positive.

- p:

  Positive whole-number exponent controlling emphasis on the shortest
  pairwise distances.

- maxit:

  Positive whole number of simulated-annealing iterations.

- temp:

  Positive initial temperature passed to
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html).

- tmax:

  Positive whole number of evaluations at each temperature.

- initial:

  Optional `n` by `q` Latin hypercube with values in `[0, 1]`.

- seed:

  Optional nonnegative whole-number seed. When supplied, the function
  restores the caller's random-number state before returning.

## Value

An object of class `magp_quantitative_design`, containing the optimized
Latin hypercube, the starting design, criterion values, minimum
distances, and search settings.

## References

Kirkpatrick, S., Gelatt, C. D., and Vecchi, M. P. (1983). Optimization
by Simulated Annealing. Science, 220, 671-680.
[doi:10.1126/science.220.4598.671](https://doi.org/10.1126/science.220.4598.671)
.

## Examples

``` r
design <- magp_quantitative_design(
  n = 8,
  q = 4,
  maxit = 500,
  seed = 1
)
design$quantity
#>      quantity_1 quantity_2 quantity_3 quantity_4
#> [1,]     0.1875     0.8125     0.5625     0.9375
#> [2,]     0.9375     0.4375     0.0625     0.6875
#> [3,]     0.5625     0.9375     0.3125     0.4375
#> [4,]     0.4375     0.1875     0.4375     0.8125
#> [5,]     0.6875     0.6875     0.9375     0.5625
#> [6,]     0.3125     0.3125     0.1875     0.0625
#> [7,]     0.8125     0.0625     0.6875     0.1875
#> [8,]     0.0625     0.5625     0.8125     0.3125
design$criterion
#> [1] 1.659959
```
