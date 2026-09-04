# Evaluate a quantitative Latin hypercube

Measures the separation of design rows under Euclidean distance. The
criterion is an inverse-distance p-norm, so smaller values indicate a
more space-filling design.

## Usage

``` r
magp_quantitative_criterion(quantity, p = 15L)
```

## Arguments

- quantity:

  Numeric matrix or data frame with values in `[0, 1]`.

- p:

  Positive whole-number exponent controlling emphasis on the shortest
  pairwise distances.

## Value

One numeric criterion value. Smaller values are preferred. A design with
duplicate rows has value `Inf`.

## Examples

``` r
quantity <- cbind(
  c(0.125, 0.375, 0.625, 0.875),
  c(0.625, 0.125, 0.875, 0.375)
)
magp_quantitative_criterion(quantity)
#> [1] 1.962421
```
