# Score candidate experiments with expected improvement

Calculates expected improvement for one or more candidate rows. Higher
values indicate candidates that offer a better combination of predicted
improvement and uncertainty.

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

  Use `"minimize"` when smaller responses are better and `"maximize"`
  when larger responses are better.

- best:

  Optional response that a new point should improve upon. When it is
  omitted, the function uses the best observed response in `object`.

- xi:

  Nonnegative improvement offset. The default is `0`. Larger values
  require a candidate to exceed `best` by more.

## Value

A nonnegative numeric vector with one expected-improvement value per row
of `newdata`.

## Details

Expected improvement uses the model's latent predictive mean and
standard error. A candidate can receive a high value because its
predicted response is good, its uncertainty is large, or both. Rows
already present in the training data receive a value of zero.
Predictions with variance below `1e-8` also receive zero.

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
