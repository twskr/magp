library(magp)

mapping_matrix_r <- function(parameters, q, mapping) {
  dimension <- if (mapping == "2d") 2L else q - 1L
  result <- matrix(0, nrow = q, ncol = dimension)
  delta <- parameters[(2L * q + 1L):length(parameters)]

  if (mapping == "2d") {
    result[2L, 2L] <- delta[1L]
    cursor <- 2L
    if (q >= 3L) {
      for (level in 3L:q) {
        for (coordinate in 1:2) {
          result[level, coordinate] <- delta[cursor]
          cursor <- cursor + 1L
        }
      }
    }
  } else {
    cursor <- 1L
    for (level in 2L:q) {
      for (coordinate in seq_len(level - 1L)) {
        result[level, coordinate] <- delta[cursor]
        cursor <- cursor + 1L
      }
    }
  }
  result
}

covariance_r <- function(X, parameters, q, tau, mapping) {
  n <- nrow(X)
  map <- mapping_matrix_r(parameters, q, mapping)
  result <- matrix(0, nrow = n, ncol = n)
  diag(result) <- tau + sum(parameters[seq_len(q)])

  for (i in seq_len(n - 1L)) {
    for (j in (i + 1L):n) {
      value <- 0
      for (component in seq_len(q)) {
        quantity_difference <- X[i, component] - X[j, component]
        first_position <- as.integer(X[i, q + component])
        second_position <- as.integer(X[j, q + component])
        mapping_difference <-
          map[first_position, ] - map[second_position, ]
        value <- value + parameters[component] * exp(
          -parameters[q + component] * quantity_difference^2 -
            sum(mapping_difference^2)
        )
      }
      result[i, j] <- value
      result[j, i] <- value
    }
  }
  result
}

time_per_call <- function(FUN, repetitions, calls_per_repetition = 1L) {
  timings <- replicate(repetitions, {
    unname(system.time({
      for (i in seq_len(calls_per_repetition)) FUN()
    })[["elapsed"]]) / calls_per_repetition
  })
  stats::median(timings)
}

benchmark_mapping <- function(mapping, X, q, tau, repetitions = 3L) {
  n_delta <- if (mapping == "2d") {
    2L * q - 3L
  } else {
    q * (q - 1L) / 2L
  }
  parameters <- c(
    rep(1, q),
    seq(0.8, 2.2, length.out = q),
    seq(-0.7, 0.7, length.out = n_delta)
  )
  mapping_code <- if (mapping == "2d") 0L else 1L

  reference <- covariance_r(X, parameters, q, tau, mapping)
  compiled <- magp:::cpp_magp_covariance(
    X,
    parameters,
    q,
    tau,
    mapping_code
  )
  difference <- max(abs(reference - compiled))
  stopifnot(difference < 1e-10)

  r_seconds <- time_per_call(
    function() covariance_r(X, parameters, q, tau, mapping),
    repetitions
  )
  rcpp_seconds <- time_per_call(
    function() magp:::cpp_magp_covariance(
      X,
      parameters,
      q,
      tau,
      mapping_code
    ),
    repetitions,
    calls_per_repetition = 100L
  )

  data.frame(
    mapping = mapping,
    rows = nrow(X),
    components = q,
    repetitions = repetitions,
    r_seconds = r_seconds,
    rcpp_seconds = rcpp_seconds,
    speedup = r_seconds / rcpp_seconds,
    max_abs_difference = difference,
    check.names = FALSE
  )
}

set.seed(20260904)
n <- 180L
q <- 8L
quantity <- matrix(stats::runif(n * q), nrow = n, ncol = q)
sequence <- t(replicate(n, sample.int(q)))
X <- cbind(quantity, sequence)
tau <- 1e-4

results <- rbind(
  benchmark_mapping("2d", X, q, tau),
  benchmark_mapping("full", X, q, tau)
)

arguments <- commandArgs(trailingOnly = TRUE)
output_directory <- if (length(arguments)) arguments[[1L]] else "."
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

write.csv(
  results,
  file.path(output_directory, "covariance-engine-results.csv"),
  row.names = FALSE
)
session_information <- sub(
  "[[:space:]]+$",
  "",
  capture.output(sessionInfo())
)
writeLines(
  session_information,
  file.path(output_directory, "covariance-engine-session-info.txt")
)
print(results)
