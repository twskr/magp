# Contributing to magp

Bug reports and focused improvements are welcome. When reporting a problem,
please include a small reproducible example, the result you expected, the
result you received, and the output of `sessionInfo()`.

For code changes:

1. Create a branch for the change.
2. Add or update tests when behavior changes.
3. Run the test suite and an R package check.
4. Open a pull request that briefly explains the change and its motivation.

You can run the tests from the package directory with:

```r
testthat::test_local()
```

Before opening a pull request, build and check the package with:

```sh
R CMD build .
R CMD check --as-cran magp_*.tar.gz
```

If an exported R or C++ interface changes, regenerate the documentation and
`Rcpp` bindings before committing the result.
