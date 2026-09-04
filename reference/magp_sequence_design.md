# Construct the sequence portion of an initial design

Uses simulated annealing to search for a sequence design with balanced
ordered adjacent pairs and well-separated rows. A neighbor is generated
by selecting one design row and exchanging two component positions, so
every candidate remains a valid permutation.

## Usage

``` r
magp_sequence_design(
  n,
  q,
  pair_weight = 0.2,
  space_weight = 0.8,
  p = 15L,
  maxit = 10000L,
  temp = 0.1,
  tmax = 10L,
  initial = NULL,
  seed = NULL
)
```

## Arguments

- n:

  Number of design runs. Must be at least two.

- q:

  Number of components. Must be at least three.

- pair_weight:

  Nonnegative weight for ordered adjacent-pair balance.

- space_weight:

  Nonnegative weight for Hamming-distance space filling.

- p:

  Positive whole-number exponent controlling emphasis on the weakest
  pair counts and the smallest distances.

- maxit:

  Positive whole number of simulated-annealing iterations.

- temp:

  Positive initial temperature passed to
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html). The default is
  scaled to the sequence-design criterion used here.

- tmax:

  Positive whole number of evaluations at each temperature.

- initial:

  Optional `n` by `q` matrix of component positions. Every row must be a
  permutation of `1:q`.

- seed:

  Optional nonnegative whole-number seed. When supplied, the function
  restores the caller's random-number state before returning.

## Value

An object of class `magp_sequence_design`, containing the optimized
sequence matrix, the starting matrix, criterion values, and search
settings. The sequence matrix is ready to use as the sequence half of a
`magp` input matrix.

## Details

This function constructs only the sequence portion of a quantitative-
sequence initial design. Use
[`magp_initial_design()`](https://twskr.github.io/magp/reference/magp_initial_design.md)
to combine it with a quantitative Latin hypercube and improve the
pairing of the two portions.

## References

Kirkpatrick, S., Gelatt, C. D., and Vecchi, M. P. (1983). Optimization
by Simulated Annealing. Science, 220, 671-680.
[doi:10.1126/science.220.4598.671](https://doi.org/10.1126/science.220.4598.671)
.

Xiao, Q., Wang, Y., Mandal, A., and Deng, X. (2024). Modeling and Active
Learning for Experiments with Quantitative-Sequence Factors. Journal of
the American Statistical Association.
[doi:10.1080/01621459.2022.2123335](https://doi.org/10.1080/01621459.2022.2123335)
.

## Examples

``` r
design <- magp_sequence_design(
  n = 8,
  q = 4,
  maxit = 500,
  seed = 1
)
design$sequence
#>      sequence_1 sequence_2 sequence_3 sequence_4
#> [1,]          4          1          2          3
#> [2,]          3          2          1          4
#> [3,]          1          4          3          2
#> [4,]          1          4          2          3
#> [5,]          3          1          4          2
#> [6,]          1          2          4          3
#> [7,]          4          2          3          1
#> [8,]          3          4          2          1
design$criterion
#> [1] 0.4509175
```
