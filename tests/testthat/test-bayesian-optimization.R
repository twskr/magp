bo_example <- function(n = 8L) {
  quantity <- rbind(
    c(0.05, 0.20, 0.80),
    c(0.20, 0.75, 0.35),
    c(0.35, 0.10, 0.55),
    c(0.50, 0.60, 0.15),
    c(0.65, 0.40, 0.95),
    c(0.80, 0.90, 0.45),
    c(0.95, 0.30, 0.65),
    c(0.42, 0.52, 0.72)
  )[seq_len(n), , drop = FALSE]
  sequence <- rbind(
    c(1, 2, 3),
    c(1, 3, 2),
    c(2, 1, 3),
    c(2, 3, 1),
    c(3, 1, 2),
    c(3, 2, 1),
    c(1, 2, 3),
    c(2, 1, 3)
  )[seq_len(n), , drop = FALSE]
  X <- cbind(quantity, sequence)
  colnames(X) <- c(
    paste0("quantity_", 1:3),
    paste0("sequence_", 1:3)
  )
  y <- -rowSums((quantity - matrix(
    c(0.2, 0.6, 0.8), nrow = n, ncol = 3, byrow = TRUE
  ))^2) - 0.01 * rowSums((sequence - matrix(
    c(1, 3, 2), nrow = n, ncol = 3, byrow = TRUE
  ))^2)
  list(X = X, y = y)
}

bo_fit <- function(mapping = "2d") {
  data <- bo_example()
  fit_function <- if (mapping == "2d") magp2d_fit else magpfull_fit
  suppressWarnings(fit_function(
    data$X,
    data$y,
    seed = 11,
    maxeval = 1
  ))
}

test_that("expected improvement matches its closed form", {
  expected <- (1 - 0) * pnorm(0.5) + 2 * dnorm(0.5)
  expect_equal(
    magp:::.magp_ei_from_prediction(1, 2, 0, "maximize", 0),
    expected,
    tolerance = 1e-14
  )
  expect_equal(
    magp:::.magp_ei_from_prediction(-1, 2, 0, "minimize", 0),
    expected,
    tolerance = 1e-14
  )
  expect_equal(
    magp:::.magp_ei_from_prediction(c(2, -1), c(0, 0), 0,
                                    "maximize", 0),
    c(0, 0)
  )
})

test_that("expected improvement uses latent model predictions", {
  fit <- bo_fit("2d")
  point <- matrix(
    c(0.12, 0.44, 0.68, 2, 3, 1),
    nrow = 1L,
    dimnames = list(NULL, colnames(fit$X))
  )
  prediction <- predict(fit, point, se.fit = TRUE, type = "latent")
  best <- max(fit$y)
  gain <- prediction$fit - best
  z <- gain / prediction$se.fit
  expected <- gain * pnorm(z) + prediction$se.fit * dnorm(z)

  expect_equal(
    magp_expected_improvement(fit, point, direction = "maximize"),
    expected,
    tolerance = 1e-12
  )
  expect_equal(
    magp_expected_improvement(
      fit, point, direction = "maximize", xi = 0.1
    ),
    magp:::.magp_ei_from_prediction(
      prediction$fit, prediction$se.fit, best, "maximize", 0.1
    ),
    tolerance = 1e-12
  )
})

test_that("expected improvement supports both fitted model classes", {
  data <- bo_example()
  point <- data$X[1:2, , drop = FALSE]
  fit_2d <- bo_fit("2d")
  fit_full <- bo_fit("full")

  expect_length(
    magp_expected_improvement(fit_2d, point, direction = "minimize"),
    2L
  )
  expect_length(
    magp_expected_improvement(fit_full, point, direction = "maximize"),
    2L
  )
  expect_equal(
    magp_expected_improvement(
      fit_2d,
      fit_2d$X[1:2, , drop = FALSE],
      direction = "maximize"
    ),
    c(0, 0)
  )
  expect_error(
    magp_expected_improvement(fit_2d, point, xi = -1),
    "xi"
  )
  expect_error(
    magp_expected_improvement(list(y = 1), point),
    "fitted"
  )
})

test_that("sequence candidates are complete or reproducibly sampled", {
  all_four <- magp:::.magp_acquisition_sequences(4L, NULL, 24L, 1L)
  expect_equal(nrow(all_four), 24L)
  expect_equal(nrow(unique(as.data.frame(all_four))), 24L)
  expect_true(all(apply(all_four, 1L, function(x) {
    identical(unname(sort(x)), 1:4)
  })))
  expect_identical(attr(all_four, "source"), "enumerated")

  first <- magp:::.magp_acquisition_sequences(6L, NULL, 10L, 8L)
  second <- magp:::.magp_acquisition_sequences(6L, NULL, 10L, 8L)
  expect_equal(first, second)
  expect_equal(nrow(first), 10L)
  expect_equal(nrow(unique(as.data.frame(first))), 10L)
  expect_identical(attr(first, "source"), "sampled")
})

test_that("next-point search is reproducible and returns an unobserved point", {
  fit <- bo_fit("2d")
  sequences <- rbind(c(1, 2, 3), c(2, 3, 1))
  first <- magp_next_point(
    fit,
    direction = "maximize",
    sequences = sequences,
    n_starts = 2,
    maxit = 3,
    seed = 20
  )
  second <- magp_next_point(
    fit,
    direction = "maximize",
    sequences = sequences,
    n_starts = 2,
    maxit = 3,
    seed = 20
  )

  expect_s3_class(first, "magp_next_point")
  expect_setequal(names(first), c(
    "call", "point", "expected_improvement", "predicted_mean",
    "predicted_standard_error", "direction", "best_observed", "xi",
    "bounds", "sequences", "sequence_source", "selected_sequence",
    "selected_start", "n_starts", "workers_requested", "workers_used",
    "execution", "diagnostics"
  ))
  expect_equal(first$point, second$point, tolerance = 0)
  expect_equal(first$expected_improvement,
               second$expected_improvement, tolerance = 0)
  expect_equal(nrow(first$diagnostics), 4L)
  expect_identical(first$sequence_source, "provided")
  expect_equal(sort(as.integer(first$point[1, 4:6])), 1:3)
  expect_false(magp:::.magp_candidate_is_observed(
    fit,
    as.matrix(first$point),
    sqrt(.Machine$double.eps)
  ))
  expect_equal(
    first$expected_improvement,
    magp_expected_improvement(
      fit,
      first$point,
      direction = "maximize"
    ),
    tolerance = 1e-12
  )
  expect_output(print(first), "expected improvement")
})

test_that("next-point controls reject invalid searches", {
  fit <- bo_fit("2d")
  expect_error(magp_next_point(fit, n_starts = 0), "n_starts")
  expect_error(magp_next_point(fit, workers = 0), "workers")
  expect_error(magp_next_point(fit, lower = -1), "prediction ranges")
  expect_error(
    magp_next_point(fit, sequences = rbind(c(1, 1, 3))),
    "permutation"
  )
  expect_error(
    magp_next_point(
      fit,
      sequences = rbind(c(1, 2, 3), c(1, 2, 3))
    ),
    "duplicated"
  )
})

test_that("parallel acquisition search uses at most two workers", {
  skip_on_cran()
  skip_if_not(
    file.exists(file.path(find.package("magp"), "Meta", "package.rds")),
    "parallel PSOCK test requires an installed package"
  )
  result <- magp_next_point(
    bo_fit("2d"),
    sequences = rbind(c(1, 2, 3)),
    n_starts = 3,
    workers = 4,
    maxit = 2,
    seed = 20
  )
  expect_equal(result$workers_requested, 4L)
  expect_equal(result$workers_used, 2L)
})

test_that("full mapping acquisition supports two components", {
  X <- rbind(
    c(0.1, 0.2, 1, 2),
    c(0.3, 0.8, 2, 1),
    c(0.5, 0.4, 1, 2),
    c(0.7, 0.6, 2, 1),
    c(0.9, 0.1, 1, 2),
    c(0.2, 0.9, 2, 1)
  )
  colnames(X) <- c(
    "quantity_1", "quantity_2", "sequence_1", "sequence_2"
  )
  y <- rowSums(X[, 1:2, drop = FALSE])
  fit <- suppressWarnings(magpfull_fit(X, y, seed = 5, maxeval = 1))
  result <- magp_next_point(
    fit,
    direction = "maximize",
    n_starts = 2,
    maxit = 2,
    seed = 6
  )

  expect_equal(nrow(result$sequences), 2L)
  expect_equal(sort(as.integer(result$point[1, 3:4])), 1:2)
  expect_true(is.finite(result$expected_improvement))
})

test_that("a supplied acquisition seed preserves the caller's random state", {
  fit <- bo_fit("2d")
  set.seed(90)
  before <- .Random.seed
  invisible(magp_next_point(
    fit,
    sequences = rbind(c(1, 2, 3)),
    n_starts = 2,
    maxit = 2,
    seed = 91
  ))
  expect_identical(.Random.seed, before)
})

test_that("acquisition bounds use the original scale of transformed inputs", {
  data <- bo_example()
  original <- data$X
  original[, 1L] <- 10 + 5 * original[, 1L]
  original[, 2L] <- -4 + 2 * original[, 2L]
  original[, 3L] <- 100 + 20 * original[, 3L]
  fit <- suppressWarnings(magp2d_fit(
    original,
    data$y,
    seed = 12,
    maxeval = 1
  ))
  result <- magp_next_point(
    fit,
    sequences = rbind(c(3, 1, 2)),
    n_starts = 2,
    maxit = 2,
    seed = 13
  )
  quantity <- as.numeric(result$point[1L, 1:3])

  expect_true(all(quantity >= fit$quantity_scaling$minimum))
  expect_true(all(quantity <= fit$quantity_scaling$maximum))
  expect_equal(result$bounds$lower, fit$quantity_scaling$minimum,
               ignore_attr = TRUE)
  expect_equal(result$bounds$upper, fit$quantity_scaling$maximum,
               ignore_attr = TRUE)
})

test_that("sequential Bayesian optimization records and refits each run", {
  data <- bo_example()
  objective <- function(quantity_1, quantity_2, quantity_3,
                        sequence_1, sequence_2, sequence_3) {
    quantity <- c(quantity_1, quantity_2, quantity_3)
    sequence <- c(sequence_1, sequence_2, sequence_3)
    -sum((quantity - c(0.2, 0.6, 0.8))^2) -
      0.01 * sum((sequence - c(1, 3, 2))^2)
  }
  result <- suppressWarnings(magp_bayes_optimize(
    objective,
    data$X,
    data$y,
    model = "2d",
    direction = "maximize",
    n_iter = 2,
    seed = 30,
    fit_control = list(maxeval = 1),
    acquisition_control = list(
      sequences = rbind(c(1, 2, 3), c(2, 3, 1)),
      n_starts = 2,
      maxit = 2
    ),
    verbose = FALSE
  ))

  expect_s3_class(result, "magp_bayes_opt")
  expect_setequal(names(result), c(
    "call", "best_point", "best_value", "best_index", "history", "model",
    "X", "y", "acquisitions", "mapping", "direction",
    "iterations_requested", "iterations_completed", "stop_reason", "xi",
    "stop_ei", "stop_patience", "fit_control", "acquisition_control"
  ))
  expect_equal(result$iterations_completed, 2L)
  expect_equal(nrow(result$X), 10L)
  expect_equal(length(result$y), 10L)
  expect_equal(nrow(result$history), 10L)
  expect_equal(length(result$acquisitions), 2L)
  expect_equal(result$model$n, 10L)
  expect_equal(result$best_value, max(result$y))
  expect_equal(length(result$best_point), 6L)
  expect_identical(result$stop_reason, "maximum iterations reached")
  expect_output(print(result), "new evaluations: 2")
})

test_that("sequential optimization accepts Score and applies early stopping", {
  data <- bo_example()
  objective <- function(quantity_1, quantity_2, quantity_3,
                        sequence_1, sequence_2, sequence_3) {
    list(Score = quantity_1 + quantity_2 + quantity_3)
  }
  result <- suppressWarnings(magp_bayes_optimize(
    objective,
    data$X,
    data$y,
    direction = "maximize",
    n_iter = 3,
    stop_ei = 1e100,
    stop_patience = 1,
    seed = 41,
    fit_control = list(maxeval = 1),
    acquisition_control = list(
      sequences = rbind(c(3, 2, 1)),
      n_starts = 2,
      maxit = 2
    ),
    verbose = FALSE
  ))

  expect_equal(result$iterations_completed, 1L)
  expect_identical(
    result$stop_reason,
    "expected improvement threshold reached"
  )
})

test_that("Bayesian optimization validates functions and control lists", {
  data <- bo_example()
  expect_error(
    magp_bayes_optimize(1, data$X, data$y, n_iter = 0),
    "FUN"
  )
  expect_error(
    magp_bayes_optimize(
      function(...) 1,
      data$X,
      data$y,
      n_iter = 0,
      fit_control = list(seed = 1)
    ),
    "unsupported"
  )
  expect_error(
    magp_bayes_optimize(
      function(...) 1,
      data$X,
      data$y,
      n_iter = 0,
      acquisition_control = list(best = 1)
    ),
    "unsupported"
  )
})

test_that("a supplied loop seed preserves the caller's random state", {
  data <- bo_example()
  set.seed(95)
  before <- .Random.seed
  result <- suppressWarnings(magp_bayes_optimize(
    function(...) 0,
    data$X,
    data$y,
    n_iter = 0,
    seed = 96,
    fit_control = list(maxeval = 1),
    verbose = FALSE
  ))
  expect_s3_class(result, "magp_bayes_opt")
  expect_identical(.Random.seed, before)
})

test_that("response columns and documented objective returns are supported", {
  data <- bo_example()
  combined <- data.frame(data$X, y = data$y, check.names = FALSE)
  result <- suppressWarnings(magp_bayes_optimize(
    FUN = function(...) 0,
    X = combined,
    n_iter = 0,
    seed = 97,
    fit_control = list(maxeval = 1),
    verbose = FALSE
  ))

  expect_equal(result$X, data$X, ignore_attr = TRUE)
  expect_equal(result$y, data$y)
  expect_true(all(result$history$Initial))
  expect_equal(nrow(result$history), nrow(data$X))

  point <- as.data.frame(data$X[1, , drop = FALSE])
  value <- magp:::.magp_bo_evaluate(
    function(..., offset) list(Value = offset),
    point,
    objective_args = list(offset = 3.5),
    context = "test"
  )
  expect_equal(value, 3.5)
  expect_error(
    magp:::.magp_bo_evaluate(
      function(...) list(result = 1),
      point,
      objective_args = list(),
      context = "test"
    ),
    "Score or Value"
  )
})
