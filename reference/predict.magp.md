# Predict outcomes from a fitted MaGP model

Returns predictions for new quantitative-sequence inputs. Plug-in
standard errors and variances can be returned with the predictions.

## Usage

``` r
# S3 method for class 'magp2d'
predict(
  object,
  newdata,
  se.fit = FALSE,
  type = c("script", "response", "latent"),
  ...
)

# S3 method for class 'magpfull'
predict(
  object,
  newdata,
  se.fit = FALSE,
  type = c("script", "response", "latent"),
  ...
)
```

## Arguments

- object:

  A fitted `magp2d` or `magpfull` object.

- newdata:

  A numeric matrix or data frame with the same quantitative and sequence
  inputs used for fitting. A response column named `y` is ignored.

- se.fit:

  Logical; if `TRUE`, return predictive standard errors and variances in
  addition to fitted values.

- type:

  Prediction convention. `"script"` uses the fitted training-row
  convention when `newdata` exactly matches a training row. `"response"`
  includes the nugget variance for a future response, and `"latent"`
  returns uncertainty for the noise-free surface.

- ...:

  Additional arguments, currently unused.

## Value

If `se.fit = FALSE`, a numeric vector of predictions. Otherwise, a list
with components `fit`, `se.fit`, `variance`, and `type`.

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
predict(fit, test[1:3, ], se.fit = TRUE, type = "response")
#> $fit
#> [1] 17.13052 29.49620 13.28978
#> 
#> $se.fit
#> [1] 30.50386 32.70545 34.35567
#> 
#> $variance
#> [1]  930.4855 1069.6463 1180.3119
#> 
#> $type
#> [1] "response"
#> 
# }
```
