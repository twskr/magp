# Construct a quantitative-sequence initial design

Builds an initial design in three steps. It first optimizes the sequence
permutations, then constructs a maximin-style Latin hypercube for the
quantitative levels, and finally aligns the two fixed designs by
permuting whole quantitative rows. The final alignment preserves both
the Latin hypercube and every sequence permutation.

## Usage

``` r
magp_initial_design(
  n,
  q,
  pair_weight = 0.2,
  sequence_space_weight = 0.8,
  quantity_weight = 0.5,
  sequence_weight = 0.5,
  p = 15L,
  sequence_maxit = 10000L,
  quantity_maxit = 10000L,
  alignment_maxit = 10000L,
  sequence_temp = 0.1,
  quantity_temp = 0.01,
  alignment_temp = 0.001,
  tmax = 10L,
  initial_sequence = NULL,
  initial_quantity = NULL,
  seed = NULL
)
```

## Arguments

- n:

  Number of design runs. Must be at least two.

- q:

  Number of components. Must be at least three.

- pair_weight:

  Nonnegative weight for ordered adjacent-pair balance in the sequence
  search.

- sequence_space_weight:

  Nonnegative weight for Hamming-distance space filling in the sequence
  search.

- quantity_weight:

  Nonnegative quantitative-distance weight in the joint criterion.

- sequence_weight:

  Nonnegative sequence-distance weight in the joint criterion.

- p:

  Positive whole-number exponent used in all three criteria.

- sequence_maxit:

  Positive whole number of sequence-search iterations.

- quantity_maxit:

  Positive whole number of quantitative-search iterations.

- alignment_maxit:

  Positive whole number of row-alignment iterations.

- sequence_temp:

  Positive initial temperature for the sequence search.

- quantity_temp:

  Positive initial temperature for the quantitative search.

- alignment_temp:

  Positive initial temperature for the alignment search.

- tmax:

  Positive whole number of evaluations at each temperature.

- initial_sequence:

  Optional `n` by `q` matrix of sequence positions.

- initial_quantity:

  Optional `n` by `q` Latin hypercube in `[0, 1]`.

- seed:

  Optional nonnegative whole-number seed. Separate deterministic seeds
  are used for the three stages, and the caller's random-number state is
  preserved.

## Value

An object of class `magp_initial_design`. Its `design` element is an `n`
by `2*q` matrix ready to use as the input to a MAGP fitting function.
The first `q` columns contain quantitative levels and the last `q`
columns contain sequence positions. The component searches and their
criterion values are retained for inspection.

## References

Xiao, Q., Wang, Y., Mandal, A., and Deng, X. (2024). Modeling and Active
Learning for Experiments with Quantitative-Sequence Factors. Journal of
the American Statistical Association.
[doi:10.1080/01621459.2022.2123335](https://doi.org/10.1080/01621459.2022.2123335)
.

## Examples

``` r
design <- magp_initial_design(
  n = 8,
  q = 4,
  sequence_maxit = 500,
  quantity_maxit = 500,
  alignment_maxit = 500,
  seed = 1
)
design$design
#>      quantity_1 quantity_2 quantity_3 quantity_4 sequence_1 sequence_2
#> [1,]     0.3125     0.8125     0.0625     0.4375          4          1
#> [2,]     0.0625     0.1875     0.4375     0.3125          3          2
#> [3,]     0.1875     0.6875     0.6875     0.8125          1          4
#> [4,]     0.4375     0.5625     0.9375     0.0625          1          4
#> [5,]     0.5625     0.0625     0.8125     0.6875          3          1
#> [6,]     0.9375     0.9375     0.5625     0.5625          1          2
#> [7,]     0.8125     0.3125     0.3125     0.1875          4          2
#> [8,]     0.6875     0.4375     0.1875     0.9375          3          4
#>      sequence_3 sequence_4
#> [1,]          2          3
#> [2,]          1          4
#> [3,]          3          2
#> [4,]          2          3
#> [5,]          4          2
#> [6,]          4          3
#> [7,]          3          1
#> [8,]          2          1
design$criteria
#>               sequence           quantitative joint_before_alignment 
#>              0.4509175              1.5794435              0.4611254 
#>                  joint 
#>              0.4521889 
```
