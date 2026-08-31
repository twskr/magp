# magp

`magp` is designed for experiments in which every component has both an amount
and a position in a sequence. It fits an additive Gaussian process that uses
both parts of the input, then returns predictions with optional uncertainty
estimates. You can choose a compact two-dimensional mapping or a full mapping
with `q - 1` dimensions. The most intensive covariance and gradient
calculations run in C++ through `Rcpp`.

## Installation

To install from a downloaded source package:

```r
install.packages(c("Rcpp", "nloptr"))
install.packages("magp_0.8.0.tar.gz", repos = NULL, type = "source")
library(magp)
```

To install from a local copy of this repository, run the following command from
its root directory:

```r
install.packages(".", repos = NULL, type = "source")
```

A C++ toolchain is required for either source installation method.

## Data format

For `q` components, the input must contain `2*q` columns:

1. the first `q` columns contain quantitative inputs;
2. the last `q` columns contain sequence positions.

Every row in the sequence columns must contain each value from `1` to `q`
exactly once. A response column named `y` may be included in the same data
frame; when it is present, the fitting functions can identify both `y` and `q`
automatically.

Quantitative columns outside `[0, 1]` are transformed by min-max scaling during
fitting. The fitted ranges are saved and used again for prediction. Inputs
already in `[0, 1]` are not changed.

## Example

```r
train <- read.table(
  system.file("extdata", "example_train.txt", package = "magp"),
  header = TRUE
)
test <- read.table(
  system.file("extdata", "example_test.txt", package = "magp"),
  header = TRUE
)

fit_2d <- magp2d_fit(train, tau = 0.001, seed = 1)
prediction_2d <- predict(fit_2d, test)
magp2d_rmse(prediction_2d, test$y)

fit_full <- magpfull_fit(train, tau = 0.001, seed = 1)
prediction_full <- predict(fit_full, test)
magp2d_rmse(prediction_full, test$y)

uncertainty <- predict(
  fit_2d,
  test[1:5, ],
  se.fit = TRUE,
  type = "response"
)
data.frame(
  prediction = uncertainty$fit,
  standard_error = uncertainty$se.fit
)
```

`tau` is a fixed nugget variance added to the covariance diagonal. It is a
variance, not a standard deviation.

## Prediction types

- `"script"` reproduces the fitted training-row convention when a new row
  exactly matches a training row.
- `"response"` includes the nugget when calculating uncertainty for a future
  response.
- `"latent"` returns uncertainty for the noise-free surface.

The reported standard errors treat the fitted covariance parameters as fixed.

## Reference

The model is described in Xiao et al. (2024),
[doi:10.1080/01621459.2022.2123335](https://doi.org/10.1080/01621459.2022.2123335).

## License

`magp` is released under the MIT License.
