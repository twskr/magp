# Fit a MaGP model with a two-dimensional sequence map

Fits an additive Gaussian process for data that combine component
amounts with a component order. Sequence positions are represented by
points in a compact two-dimensional latent map.

## Usage

``` r
magp2d_fit(
  X,
  y = NULL,
  q = NULL,
  tau = 0.001,
  maxeval = 500,
  xtol_rel = 1e-05,
  lb_sigma = 10,
  ub_sigma = 1000,
  lb_theta = 0.5,
  ub_theta = 1000,
  lb_delta = -1,
  ub_delta = 1,
  seed = NULL,
  n_starts = 1,
  workers = 1
)
```

## Arguments

- X:

  A numeric matrix or data frame. The first `q` input columns contain
  quantitative levels and the next `q` columns contain sequence
  positions. Each sequence row must be a permutation of `1:q`. A column
  named `y` may be included as the response.

- y:

  An optional numeric response vector. It may be omitted when `X`
  contains a response column named `y`.

- q:

  The number of components. It is inferred from the number of input
  columns when omitted. The two-dimensional model requires at least
  three components; the full model requires at least two.

- tau:

  A fixed nonnegative nugget variance added to the covariance diagonal.

- maxeval:

  Maximum number of objective evaluations used by `nloptr`.

- xtol_rel:

  Relative parameter tolerance used by `nloptr`.

- lb_sigma, ub_sigma:

  Lower and upper bounds for the additive variance parameters.

- lb_theta, ub_theta:

  Lower and upper bounds for the quantitative correlation parameters.

- lb_delta, ub_delta:

  Lower and upper bounds for the mapping parameters.

- seed:

  An optional nonnegative integer used to generate the initial parameter
  vectors. The first start retains the result produced by this seed when
  `n_starts = 1`.

- n_starts:

  Number of independent parameter starts. The fitted object contains the
  result with the lowest objective among the converged starts.

- workers:

  Number of local worker processes. Values greater than one use a socket
  cluster and are capped at `n_starts`.

## Value

An object of class `magp2d`.

## Details

The covariance is a sum of component-specific terms. Each term combines
the distance between two quantitative levels with the distance between
their mapped sequence positions. The variance, correlation, and mapping
parameters are estimated together with bounded optimization.

Quantitative columns outside `[0, 1]` are transformed with column-wise
min-max scaling. Their training ranges are stored in the fitted object
and reused for prediction. Columns already in `[0, 1]` are left
unchanged.

When several starts are requested, each start receives a separate seed.
Supplying `seed` makes the starts and the selected fit reproducible for
both sequential and parallel execution. Start-level objective values,
convergence codes, warnings, and errors are stored in
`fit$multistart$starts`.

## Examples

``` r
# \donttest{
train <- read.table(
  system.file("extdata", "example_train.txt", package = "magp"),
  header = TRUE
)
fit <- magp2d_fit(train, seed = 1, n_starts = 2)
fit
#> <magp2d model>
#>   components: 4 
#>   training runs: 31 
#>   nugget: 0.001 
#>   fitted mean: 15.0858 
#>   objective: 205.1974 
#>   optimizer: converged (status 4 )
#>   starts: 2 (best 2 )
#>   execution: sequential with 1 worker(s)
# }
```
