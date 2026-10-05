workflow_objective <- function(quantity_1, quantity_2, quantity_3,
                               sequence_1, sequence_2, sequence_3) {
  quantity <- c(quantity_1, quantity_2, quantity_3)
  sequence <- c(sequence_1, sequence_2, sequence_3)
  -sum((quantity - c(0.2, 0.6, 0.8))^2) -
    0.01 * sum((sequence - c(1, 3, 2))^2)
}

workflow_design_control <- list(
  sequence_maxit = 30,
  quantity_maxit = 30,
  alignment_maxit = 30
)

workflow_latin_hypercube_is_valid <- function(quantity) {
  n <- nrow(quantity)
  strata <- floor(pmin(quantity, 1 - .Machine$double.eps) * n)
  all(apply(strata, 2L, function(column) {
    identical(sort(as.integer(column)), 0:(n - 1L))
  }))
}

test_that("from-scratch optimization connects every workflow stage", {
  output <- capture.output({
    result <- suppressWarnings(magp_bayes_optimize_from_scratch(
      workflow_objective,
      n_initial = 6,
      q = 3,
      direction = "maximize",
      n_iter = 2,
      seed = 101,
      design_control = workflow_design_control,
      fit_control = list(maxeval = 2),
      acquisition_control = list(n_starts = 2, maxit = 2),
      verbose = TRUE
    ))
  })

  expected_initial <- apply(
    result$initial_design$design,
    1L,
    function(row) do.call(workflow_objective, as.list(row))
  )
  expect_s3_class(result, "magp_bayes_opt")
  expect_s3_class(result$initial_design, "magp_initial_design")
  expect_setequal(names(result), c(
    "call", "best_point", "best_value", "best_index", "history", "model",
    "X", "y", "acquisitions", "mapping", "direction",
    "iterations_requested", "iterations_completed", "stop_reason", "xi",
    "stop_ei", "stop_patience", "fit_control", "acquisition_control",
    "initial_design", "initial_response", "initial_evaluations",
    "design_control", "started_from_initial_design"
  ))
  expect_true(result$started_from_initial_design)
  expect_equal(result$initial_evaluations, 6L)
  expect_equal(result$initial_response, expected_initial, tolerance = 1e-14)
  expect_equal(result$iterations_completed, 2L)
  expect_equal(result$model$n, 8L)
  expect_equal(nrow(result$history), 8L)
  expect_equal(sum(result$history$Initial), 6L)
  expect_equal(length(result$acquisitions), 2L)
  expect_false(anyDuplicated(as.data.frame(result$X)) > 0L)
  expect_true(workflow_latin_hypercube_is_valid(
    result$initial_design$quantity
  ))
  expect_true(all(apply(result$initial_design$sequence, 1L, function(row) {
    identical(sort(as.integer(row)), 1:3)
  })))
  expect_true(any(grepl("initial run 1", output, fixed = TRUE)))
  expect_true(any(grepl("iteration 1", output, fixed = TRUE)))
  expect_output(print(result), "initial design runs: 6")
})

test_that("from-scratch optimization supports full mapping and fixed args", {
  objective <- function(quantity_1, quantity_2, quantity_3,
                        sequence_1, sequence_2, sequence_3, offset) {
    list(Score = quantity_1 + quantity_2 + quantity_3 + offset)
  }
  result <- suppressWarnings(magp_bayes_optimize_from_scratch(
    objective,
    n_initial = 6,
    q = 3,
    model = "full",
    direction = "minimize",
    n_iter = 0,
    seed = 102,
    design_control = workflow_design_control,
    fit_control = list(maxeval = 1),
    objective_args = list(offset = 2),
    verbose = FALSE
  ))

  expected <- rowSums(result$initial_design$quantity) + 2
  expect_s3_class(result$model, "magpfull")
  expect_equal(result$initial_response, expected, tolerance = 1e-14)
  expect_equal(result$y, expected, tolerance = 1e-14)
  expect_equal(result$best_value, min(expected))
  expect_equal(result$iterations_completed, 0L)
})

test_that("from-scratch optimization is reproducible and preserves RNG", {
  stochastic_objective <- function(quantity_1, quantity_2, quantity_3,
                                   sequence_1, sequence_2, sequence_3) {
    quantity_1 + stats::runif(1)
  }
  set.seed(103)
  state_before <- .Random.seed
  first <- suppressWarnings(magp_bayes_optimize_from_scratch(
    stochastic_objective,
    n_initial = 5,
    q = 3,
    n_iter = 0,
    seed = 104,
    design_control = workflow_design_control,
    fit_control = list(maxeval = 1),
    verbose = FALSE
  ))
  state_after_first <- .Random.seed
  second <- suppressWarnings(magp_bayes_optimize_from_scratch(
    stochastic_objective,
    n_initial = 5,
    q = 3,
    n_iter = 0,
    seed = 104,
    design_control = workflow_design_control,
    fit_control = list(maxeval = 1),
    verbose = FALSE
  ))

  expect_identical(state_after_first, state_before)
  expect_identical(.Random.seed, state_before)
  expect_equal(first$initial_design$design, second$initial_design$design)
  expect_equal(first$initial_response, second$initial_response)
  expect_equal(first$best_point, second$best_point)
  expect_equal(first$best_value, second$best_value)
})

test_that("from-scratch controls fail before expensive work begins", {
  expect_error(
    magp_bayes_optimize_from_scratch(1, 6, 3),
    "FUN"
  )
  expect_error(
    magp_bayes_optimize_from_scratch(workflow_objective, 6, 2),
    "q"
  )
  expect_error(
    magp_bayes_optimize_from_scratch(
      workflow_objective, 6, 3,
      design_control = list(seed = 1)
    ),
    "unsupported"
  )
  expect_error(
    magp_bayes_optimize_from_scratch(
      workflow_objective, 6, 3,
      fit_control = list(seed = 1)
    ),
    "unsupported"
  )
  expect_error(
    magp_bayes_optimize_from_scratch(
      workflow_objective, 6, 3,
      acquisition_control = list(best = 1)
    ),
    "unsupported"
  )
  expect_error(
    magp_bayes_optimize_from_scratch(
      workflow_objective, 6, 3,
      objective_args = list(1)
    ),
    "nonempty names"
  )
  expect_error(
    magp_bayes_optimize_from_scratch(
      workflow_objective, 6, 3,
      objective_args = list(quantity_1 = 1)
    ),
    "duplicates input"
  )
})

test_that("initial objective errors identify the failed run", {
  broken <- function(...) stop("objective unavailable")
  expect_error(
    magp_bayes_optimize_from_scratch(
      broken,
      n_initial = 4,
      q = 3,
      n_iter = 0,
      seed = 105,
      design_control = list(
        sequence_maxit = 1,
        quantity_maxit = 1,
        alignment_maxit = 1
      ),
      fit_control = list(maxeval = 1),
      verbose = FALSE
    ),
    "initial run 1: objective unavailable"
  )
})
