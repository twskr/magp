# Contributing to magp

Bug reports and focused improvements are welcome. When reporting a
problem, please include a small reproducible example, the result you
expected, the result you received, and the output of
[`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html).

For code changes:

1.  Create a branch for the change.
2.  Add or update tests when behavior changes.
3.  Run the test suite and an R package check.
4.  Open a pull request that briefly explains the change and its
    motivation.

The full test suite includes checks that start fresh R worker processes.
Run it through an installed-package check from the package directory:

``` r

devtools::check()
```

You can also build and check it directly:

``` sh
R CMD build .
R CMD check --as-cran magp_*.tar.gz
```

If an exported R or C++ interface changes, regenerate the documentation
and `Rcpp` bindings before committing the result.
