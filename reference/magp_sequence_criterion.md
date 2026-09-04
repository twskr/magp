# Evaluate a sequence initial design

Measures two properties of a sequence design: how evenly ordered
adjacent component pairs are represented and how well separated the
design rows are under Hamming distance. Smaller values indicate a better
design.

## Usage

``` r
magp_sequence_criterion(
  sequence,
  pair_weight = 0.2,
  space_weight = 0.8,
  p = 15L
)
```

## Arguments

- sequence:

  Numeric matrix or data frame. Every row must be a permutation of
  `1:q`.

- pair_weight:

  Nonnegative weight for ordered adjacent-pair balance.

- space_weight:

  Nonnegative weight for Hamming-distance space filling.

- p:

  Positive whole-number exponent controlling emphasis on the weakest
  pair counts and the smallest distances.

## Value

One numeric criterion value. Smaller values are preferred.

## Details

Each row uses the same sequence format as the fitting functions. The
value in column `j` is the position assigned to component `j`, and every
row must be a permutation of `1:q`.

## References

Xiao, Q., Wang, Y., Mandal, A., and Deng, X. (2024). Modeling and Active
Learning for Experiments with Quantitative-Sequence Factors. Journal of
the American Statistical Association.
[doi:10.1080/01621459.2022.2123335](https://doi.org/10.1080/01621459.2022.2123335)
.

## Examples

``` r
sequence <- rbind(
  c(1, 2, 3, 4),
  c(3, 1, 4, 2),
  c(2, 4, 1, 3),
  c(4, 3, 2, 1)
)
magp_sequence_criterion(sequence)
#> [1] 0.5300508
```
