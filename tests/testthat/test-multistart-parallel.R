multistart_example <- function() {
  train <- read.table(
    system.file("extdata", "example_train.txt", package = "magp"),
    header = TRUE
  )
  columns <- c("A", "B", "C", "D", "a", "b", "c", "d")
  list(
    X = as.matrix(train[seq_len(8L), columns]),
    y = train$y[seq_len(8L)]
  )
}

test_that("multi-start controls are validated", {
  data <- multistart_example()

  expect_error(
    magp2d_fit(data$X, data$y, n_starts = 0),
    "n_starts"
  )
  expect_error(
    magp2d_fit(data$X, data$y, n_starts = 1.5),
    "n_starts"
  )
  expect_error(
    magp2d_fit(data$X, data$y, workers = 0),
    "workers"
  )
  expect_error(
    magpfull_fit(data$X, data$y, workers = NA_real_),
    "workers"
  )
})

test_that("one start retains the established seeded fit", {
  data <- multistart_example()
  fit <- suppressWarnings(magp2d_fit(
    data$X, data$y,
    seed = 17,
    maxeval = 1
  ))
  direct <- suppressWarnings(magp:::.magp_fast_fit(
    data$X, data$y, NULL,
    0.001, 1, 1e-5,
    10, 1000, 0.5, 1000,
    -1, 1, 17, "2d", quote(direct_fit())
  ))

  expect_equal(fit$solution, direct$solution, tolerance = 0)
  expect_equal(fit$objective, direct$objective, tolerance = 0)
  expect_equal(fit$multistart$n_starts, 1L)
  expect_equal(fit$multistart$workers_used, 1L)
  expect_equal(fit$multistart$best_start, 1L)
  expect_equal(fit$multistart$starts$seed, 17L)
})

test_that("seeded sequential starts are reproducible", {
  data <- multistart_example()
  first <- suppressWarnings(magp2d_fit(
    data$X, data$y,
    seed = 23,
    n_starts = 3,
    workers = 1,
    maxeval = 1
  ))
  second <- suppressWarnings(magp2d_fit(
    data$X, data$y,
    seed = 23,
    n_starts = 3,
    workers = 1,
    maxeval = 1
  ))

  expect_equal(first$solution, second$solution, tolerance = 0)
  expect_equal(first$objective, second$objective, tolerance = 0)
  expect_equal(first$multistart$starts, second$multistart$starts)
  expect_equal(first$multistart$n_starts, 3L)
  expect_equal(first$multistart$workers_used, 1L)
  expect_identical(first$multistart$mode, "sequential")
  expect_identical(as.character(first$call[[1L]]), "magp2d_fit")
  expect_output(print(first), "starts: 3")
  expect_output(print(first), "execution: sequential")
  expect_equal(
    first$objective,
    min(first$multistart$starts$objective, na.rm = TRUE),
    tolerance = 0
  )
})

test_that("full mapping supports multiple starts", {
  data <- multistart_example()
  fit <- suppressWarnings(magpfull_fit(
    data$X, data$y,
    seed = 31,
    n_starts = 2,
    workers = 1,
    maxeval = 1
  ))

  expect_s3_class(fit, "magpfull")
  expect_equal(fit$multistart$n_starts, 2L)
  expect_equal(nrow(fit$multistart$starts), 2L)
  expect_true(all(is.finite(fit$multistart$starts$objective)))
  expect_true(fit$multistart$best_start %in% 1:2)
})

test_that("the lowest converged objective is selected", {
  data <- multistart_example()
  fit <- magp2d_fit(
    data$X, data$y,
    seed = 9,
    n_starts = 3,
    workers = 1,
    maxeval = 500
  )
  starts <- fit$multistart$starts
  converged <- which(starts$converged)
  expected_best <- converged[which.min(starts$objective[converged])]

  expect_true(length(converged) > 0L)
  expect_equal(
    fit$objective,
    min(starts$objective[converged]),
    tolerance = 1e-12
  )
  expect_equal(fit$multistart$best_start, expected_best)
})

test_that("parallel and sequential starts agree for both models", {
  skip_on_cran()
  skip_if_not(
    file.exists(file.path(find.package("magp"), "Meta", "package.rds")),
    "parallel PSOCK test requires an installed package"
  )
  data <- multistart_example()
  sequential <- suppressWarnings(magp2d_fit(
    data$X, data$y,
    seed = 47,
    n_starts = 3,
    workers = 1,
    maxeval = 1
  ))
  parallel_fit <- suppressWarnings(magp2d_fit(
    data$X, data$y,
    seed = 47,
    n_starts = 3,
    workers = 4,
    maxeval = 1
  ))

  expect_equal(parallel_fit$solution, sequential$solution, tolerance = 0)
  expect_equal(parallel_fit$objective, sequential$objective, tolerance = 0)
  expect_equal(
    parallel_fit$multistart$starts,
    sequential$multistart$starts
  )
  expect_identical(parallel_fit$multistart$mode, "PSOCK")
  expect_equal(parallel_fit$multistart$workers_used, 2L)

  sequential_full <- suppressWarnings(magpfull_fit(
    data$X, data$y,
    seed = 53,
    n_starts = 2,
    workers = 1,
    maxeval = 1
  ))
  parallel_full <- suppressWarnings(magpfull_fit(
    data$X, data$y,
    seed = 53,
    n_starts = 2,
    workers = 4,
    maxeval = 1
  ))

  expect_equal(
    parallel_full$solution,
    sequential_full$solution,
    tolerance = 0
  )
  expect_equal(
    parallel_full$objective,
    sequential_full$objective,
    tolerance = 0
  )
  expect_equal(
    parallel_full$multistart$starts,
    sequential_full$multistart$starts
  )
  expect_equal(parallel_full$multistart$workers_used, 2L)
})
