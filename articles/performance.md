# Computational benchmark

The covariance calculation is the most frequently repeated low-level
operation during fitting and prediction. `magp` implements it in C++
through `Rcpp`. This page records a direct comparison with an
independent R implementation of the same calculation.

## Reproduce the benchmark

Clone the repository, install the development package, and run:

``` sh
Rscript inst/benchmarks/benchmark-covariance-engine.R benchmark-output
```

The script creates fixed synthetic inputs, times both implementations,
checks their largest absolute numerical difference, and records
[`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html). The
complete script is included in the repository so that the comparison can
be repeated on other computers.

## Recorded result

``` r

results <- read.csv(system.file(
  "benchmarks",
  "covariance-engine-results.csv",
  package = "magp"
))
knitr::kable(
  results,
  digits = 4,
  caption = "Median elapsed time per covariance calculation."
)
```

| mapping | rows | components | repetitions | r_seconds | rcpp_seconds | speedup | max_abs_difference |
|:---|---:|---:|---:|---:|---:|---:|---:|
| 2d | 180 | 8 | 3 | 0.125 | 0.0015 | 85.6164 | 0 |
| full | 180 | 8 | 3 | 0.151 | 0.0018 | 83.4254 | 0 |

Median elapsed time per covariance calculation. {.table
style="width:100%;"}

`speedup` is the R time divided by the Rcpp time. Timing depends on the
processor, R version, compiler, and system load, so the recorded values
should not be treated as a universal performance guarantee. The
`max_abs_difference` column is the direct numerical agreement check for
the same inputs and parameters.

The benchmark session details are stored beside the result file:

``` r

session_file <- system.file(
  "benchmarks",
  "covariance-engine-session-info.txt",
  package = "magp"
)
cat(paste(readLines(session_file), collapse = "\n"))
#> R version 4.5.2 (2025-10-31)
#> Platform: aarch64-apple-darwin20
#> Running under: macOS Tahoe 26.5.2
#> 
#> Matrix products: default
#> BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
#> 
#> locale:
#> [1] C.UTF-8/C.UTF-8/C.UTF-8/C/C.UTF-8/C.UTF-8
#> 
#> time zone: America/New_York
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base
#> 
#> other attached packages:
#> [1] magp_0.12.0.9000
#> 
#> loaded via a namespace (and not attached):
#> [1] compiler_4.5.2 tools_4.5.2    Rcpp_1.1.2
```
