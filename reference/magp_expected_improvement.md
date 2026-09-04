# Calculate expected improvement for a fitted MaGP model

Evaluates the expected improvement acquisition function at one or more
quantitative-sequence inputs. Predictive uncertainty is taken from the
latent, noise-free surface.

## Usage

``` r
magp_expected_improvement(
  object,
  newdata,
  direction = c("minimize", "maximize"),
  best = NULL,
  xi = 0
)
```

## Arguments

- object:

  A fitted `magp2d` or `magpfull` model.

- newdata:

  A numeric matrix or data frame accepted by the model's
  [`stats::predict()`](https://rdrr.io/r/stats/predict.html) method.

- direction:

  Whether the objective is being minimized or maximized.

- best:

  Optional finite reference value. When omitted, the smallest or largest
  observed response is used according to `direction`.

- xi:

  A nonnegative exploration offset. Larger values require a greater
  improvement over `best` before favoring exploitation.

## Value

A nonnegative numeric vector with one value per row of `newdata`.

## Details

Let `m` and `s` denote the latent predictive mean and standard error.
For maximization, the improvement is `m - best - xi`; for minimization,
it is `best - m - xi`. If `s` is positive, expected improvement is
`improvement * pnorm(z) + s * dnorm(z)`, where `z = improvement / s`.
Values with predictive variance below `1e-8`, and inputs already present
in the training data, receive expected improvement zero.

## References

Jones, D. R., Schonlau, M., and Welch, W. J. (1998). Efficient Global
Optimization of Expensive Black-Box Functions. Journal of Global
Optimization, 13, 455-492.
[doi:10.1023/A:1008306431147](https://doi.org/10.1023/A%3A1008306431147)
.

## Examples

``` r
# \donttest{
train <- read.table(
  system.file("extdata", "example_train.txt", package = "magp"),
  header = TRUE
)
test <- read.table(
  system.file("extdata", "example_test.txt", package = "magp"),
  header = TRUE
)
fit <- magp2d_fit(train, seed = 1)
magp_expected_improvement(
  fit,
  test[1:3, ],
  direction = "minimize"
)
#> [1] 1.0641766 0.5910628 2.1282648
# }
```
