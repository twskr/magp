# Evaluate a complete quantitative-sequence initial design

Combines Euclidean distance for the quantitative portion with Hamming
distance for the sequence portion. Smaller values indicate better joint
separation of the design runs.

## Usage

``` r
magp_joint_criterion(
  quantity,
  sequence,
  quantity_weight = 0.5,
  sequence_weight = 0.5,
  p = 15L
)
```

## Arguments

- quantity:

  Numeric matrix or data frame with values in `[0, 1]`.

- sequence:

  Numeric matrix or data frame with the same dimensions as `quantity`.
  Every row must be a permutation of `1:q`.

- quantity_weight:

  Nonnegative weight for quantitative distance.

- sequence_weight:

  Nonnegative weight for sequence distance.

- p:

  Positive whole-number exponent controlling emphasis on the least
  separated pairs of runs.

## Value

One numeric criterion value. Smaller values are preferred.

## References

Xiao, Q., Wang, Y., Mandal, A., and Deng, X. (2024). Modeling and Active
Learning for Experiments with Quantitative-Sequence Factors. Journal of
the American Statistical Association.
[doi:10.1080/01621459.2022.2123335](https://doi.org/10.1080/01621459.2022.2123335)
.

## Examples

``` r
quantity <- cbind(
  c(0.125, 0.375, 0.625, 0.875),
  c(0.625, 0.125, 0.875, 0.375),
  c(0.375, 0.875, 0.125, 0.625),
  c(0.875, 0.625, 0.375, 0.125)
)
sequence <- rbind(
  c(1, 2, 3, 4),
  c(3, 1, 4, 2),
  c(2, 4, 1, 3),
  c(4, 3, 2, 1)
)
magp_joint_criterion(quantity, sequence)
#> [1] 0.3278273
```
