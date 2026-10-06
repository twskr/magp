# Independent sequence-criterion calculation in execution-order form.
# Non-self-inverse permutations make a missing positions-to-orders conversion
# visible in the adjacent-pair term.
sfta_reference_criterion <- function(sequence, pair_weight = 0.2,
                                     space_weight = 0.8, p = 15L) {
  orders <- t(apply(sequence, 1L, order))
  n <- nrow(orders)
  q <- ncol(orders)
  counts <- matrix(0, q, q)
  for (i in seq_len(n)) {
    for (j in seq_len(q - 1L)) {
      counts[orders[i, j], orders[i, j + 1L]] <-
        counts[orders[i, j], orders[i, j + 1L]] + 1
    }
  }
  pair_term <- sum((1 + counts[row(counts) != col(counts)])^(-p))
  space_term <- 0
  for (i in seq_len(n - 1L)) {
    for (j in seq.int(i + 1L, n)) {
      distance <- sum(orders[i, ] != orders[j, ])
      space_term <- space_term + (1 + distance)^(-p)
    }
  }
  ((pair_weight * pair_term + space_weight * space_term) /
     (pair_weight + space_weight))^(1 / p)
}

sfta_test_control <- function() {
  list(nstarts = 2L, ncalibrate = 20L, nrounds = 4L,
       max_proposals = 100L)
}

sfta_valid_positions <- function(sequence) {
  all(apply(sequence, 1L, function(row) {
    identical(sort(as.integer(row)), seq_len(ncol(sequence)))
  }))
}

sfta_valid_latin <- function(quantity) {
  n <- nrow(quantity)
  all(apply(floor(pmin(quantity, 1 - .Machine$double.eps) * n),
            2L, function(column) {
    identical(sort(as.integer(column)), seq.int(0L, n - 1L))
  }))
}

test_that("random sequence designs retain their generated or supplied start", {
  generated <- magp_sequence_design(
    n = 8, q = 4, method = "random", seed = 611
  )
  repeated <- magp_sequence_design(
    n = 8, q = 4, method = "random", seed = 611
  )
  expect_s3_class(generated, "magp_sequence_design")
  expect_identical(generated$method_key, "random")
  expect_identical(generated$method, "random sampling")
  expect_equal(generated$sequence, repeated$sequence)
  expect_equal(generated$sequence, generated$initial_sequence)
  expect_equal(generated$criterion, generated$initial_criterion)
  expect_equal(generated$improvement, 0)
  expect_true(sfta_valid_positions(generated$sequence))
  expect_equal(generated$criterion,
               sfta_reference_criterion(generated$sequence),
               tolerance = 1e-13)

  initial <- rbind(c(2, 3, 1, 4), c(3, 4, 2, 1), c(1, 4, 3, 2))
  supplied <- magp_sequence_design(
    n = 3, q = 4, initial = initial, method = "random", seed = 612
  )
  expect_equal(unname(supplied$sequence), initial)
  expect_equal(supplied$sequence, supplied$initial_sequence)
  expect_equal(supplied$improvement, 0)
})

test_that("SFTA preserves permutations, the baseline, and reproducibility", {
  initial <- rbind(
    c(2, 3, 1, 4), c(3, 4, 2, 1), c(1, 4, 3, 2),
    c(4, 2, 1, 3), c(2, 1, 4, 3), c(3, 2, 4, 1)
  )
  first <- magp_sequence_design(
    n = 6, q = 4, initial = initial, method = "sfta", maxit = 40,
    sfta_control = sfta_test_control(), seed = 613
  )
  second <- magp_sequence_design(
    n = 6, q = 4, initial = initial, method = "sfta", maxit = 40,
    sfta_control = sfta_test_control(), seed = 613
  )
  expect_identical(first$method_key, "sfta")
  expect_identical(first$method, "space-filling threshold accepting")
  expect_equal(first$sequence, second$sequence)
  expect_equal(first$criterion, second$criterion)
  expect_equal(first$sfta, second$sfta)
  expect_equal(unname(first$initial_sequence), initial)
  expect_equal(dim(first$sequence), c(6L, 4L))
  expect_identical(colnames(first$sequence), paste0("sequence_", 1:4))
  expect_true(sfta_valid_positions(first$sequence))
  expect_lte(first$criterion, first$initial_criterion + 1e-13)
  expect_equal(first$initial_criterion, sfta_reference_criterion(initial),
               tolerance = 1e-13)
  expect_equal(first$criterion, sfta_reference_criterion(first$sequence),
               tolerance = 1e-13)
  expect_equal(first$improvement,
               first$initial_criterion - first$criterion,
               tolerance = 1e-13)

  expect_type(first$sfta, "list")
  expect_true(all(c("controls", "thresholds", "phase1_criterion",
                    "phase1_attempts", "phase1_fallbacks", "accepted",
                    "evaluations") %in% names(first$sfta)))
  expect_true(all(is.finite(first$sfta$thresholds)))
  expect_true(all(first$sfta$thresholds >= 0))
  expect_true(all(is.finite(first$sfta$phase1_criterion)))
  expect_lte(first$criterion, min(first$sfta$phase1_criterion) + 1e-13)
  expect_true(all(first$sfta$accepted >= 0))
  expect_lte(sum(first$sfta$accepted), 40L)
})

test_that("SFTA criterion agrees with independent pair and spacing terms", {
  initial <- rbind(
    c(2, 3, 1, 4), c(3, 4, 2, 1), c(1, 4, 3, 2),
    c(4, 2, 1, 3), c(2, 1, 4, 3), c(3, 2, 4, 1)
  )
  for (p in c(1L, 15L)) {
    for (weights in list(c(1, 0), c(0, 1), c(2, 3))) {
      fit <- magp_sequence_design(
        n = 6, q = 4, initial = initial, method = "sfta", maxit = 20,
        pair_weight = weights[1L], space_weight = weights[2L], p = p,
        sfta_control = list(nstarts = 2L, ncalibrate = 6L,
                            nrounds = 2L, max_proposals = 50L),
        seed = 614
      )
      expect_true(sfta_valid_positions(fit$sequence))
      expect_equal(fit$criterion,
                   sfta_reference_criterion(fit$sequence,
                                            weights[1L], weights[2L], p),
                   tolerance = 1e-12)
      expect_lte(fit$criterion, fit$initial_criterion + 1e-12)
    }
  }
})

test_that("space-filling sampling terminates after exhausting permutations", {
  for (n in c(6L, 8L)) {
    fit <- magp_sequence_design(
      n = n, q = 3, method = "sfta", maxit = 8,
      sfta_control = list(nstarts = 2L, ncalibrate = 4L,
                          nrounds = 2L, max_proposals = 1L),
      seed = 615
    )
    expect_equal(dim(fit$sequence), c(n, 3L))
    expect_true(sfta_valid_positions(fit$sequence))
    expect_true(is.finite(fit$criterion))
    expect_true(any(fit$sfta$phase1_fallbacks > 0))
    expect_equal(fit$criterion, sfta_reference_criterion(fit$sequence),
                 tolerance = 1e-13)
  }
})

test_that("new sequence methods preserve the caller's random state", {
  for (method in c("random", "sfta")) {
    set.seed(616)
    before <- .Random.seed
    invisible(magp_sequence_design(
      n = 6, q = 3, maxit = 40, method = method, seed = 617,
      sfta_control = if (method == "sfta") sfta_test_control() else list()
    ))
    expect_identical(.Random.seed, before)
  }
})

test_that("legacy default remains the same as explicitly selected SANN", {
  default <- magp_sequence_design(n = 6, q = 4, maxit = 40, seed = 618)
  explicit <- magp_sequence_design(
    n = 6, q = 4, maxit = 40, seed = 618, method = "sann"
  )
  expect_identical(default$method_key, "sann")
  expect_identical(default$method, "simulated annealing")
  expect_equal(default$sequence, explicit$sequence)
  expect_equal(default$initial_sequence, explicit$initial_sequence)
  expect_equal(default$criterion, explicit$criterion)

  args <- list(n = 6, q = 4, sequence_maxit = 40, quantity_maxit = 20,
               alignment_maxit = 20, seed = 619)
  default_design <- do.call(magp_initial_design, args)
  explicit_design <- do.call(magp_initial_design,
                             c(args, list(sequence_method = "sann")))
  expect_equal(default_design$design, explicit_design$design)
  expect_equal(default_design$criteria, explicit_design$criteria)
})

test_that("sequence selection preserves the quantitative design and layout", {
  designs <- lapply(c("sann", "random", "sfta"), function(method) {
    magp_initial_design(
      n = 8, q = 4, sequence_maxit = 40, quantity_maxit = 30,
      alignment_maxit = 30, seed = 620, sequence_method = method,
      sfta_control = if (method == "sfta") sfta_test_control() else list()
    )
  })
  for (i in seq_along(designs)) {
    design <- designs[[i]]
    expect_s3_class(design, "magp_initial_design")
    expect_equal(design$quantity_search$quantity,
                 designs[[1L]]$quantity_search$quantity)
    expect_identical(design$sequence_search$method_key,
                      c("sann", "random", "sfta")[[i]])
    expect_equal(design$sequence, design$sequence_search$sequence)
    expect_true(sfta_valid_latin(design$quantity))
    expect_true(sfta_valid_positions(design$sequence))
    expect_equal(design$quantity,
                 design$quantity_search$quantity[design$alignment$row_order,
                                                 , drop = FALSE])
    expect_equal(design$design, cbind(design$quantity, design$sequence))
    expect_identical(colnames(design$design),
                      c(paste0("quantity_", 1:4), paste0("sequence_", 1:4)))
    expect_equal(design$criteria[["sequence"]],
                 sfta_reference_criterion(design$sequence), tolerance = 1e-13)
  }
})

test_that("SFTA controls and method selection are validated", {
  expect_error(magp_sequence_design(n = 4, q = 3, method = "unknown"))
  expect_error(magp_initial_design(n = 4, q = 3,
                                   sequence_method = "unknown"))
  invalid_controls <- list(
    list(unknown = 1), list(1), list(nstarts = 1, nstarts = 2),
    list(nstarts = 0), list(nstarts = 1.5), list(ncalibrate = NA_real_),
    list(nrounds = 41L), list(max_proposals = 0),
    list(nstarts = as.double(.Machine$integer.max) + 1)
  )
  for (control in invalid_controls) {
    expect_error(magp_sequence_design(
      n = 4, q = 3, method = "sfta", maxit = 40,
      sfta_control = control
    ))
  }
  for (method in c("sann", "random")) {
    expect_error(magp_sequence_design(
      n = 4, q = 3, method = method, sfta_control = list(nstarts = 2L)
    ))
  }
})

test_that("from-scratch design controls propagate to initial sequence search", {
  objective <- function(quantity_1, quantity_2, quantity_3,
                        sequence_1, sequence_2, sequence_3) {
    sum((c(quantity_1, quantity_2, quantity_3) - c(0.2, 0.6, 0.8))^2) +
      0.1 * sum((c(sequence_1, sequence_2, sequence_3) - c(1, 3, 2))^2)
  }
  for (method in c("random", "sfta")) {
    controls <- list(
      sequence_method = method, sequence_maxit = 40,
      quantity_maxit = 20, alignment_maxit = 20,
      sfta_control = if (method == "sfta") sfta_test_control() else list()
    )
    result <- suppressWarnings(magp_bayes_optimize_from_scratch(
      FUN = objective, n_initial = 5, q = 3, n_iter = 0,
      design_control = controls, fit_control = list(maxeval = 2),
      seed = 621, verbose = FALSE
    ))
    expected_y <- apply(result$initial_design$design, 1L, function(row) {
      do.call(objective, as.list(row))
    })
    expect_s3_class(result, "magp_bayes_opt")
    expect_identical(result$initial_design$sequence_search$method_key, method)
    expect_equal(result$initial_response, expected_y, tolerance = 1e-13)
    expect_equal(result$initial_evaluations, 5L)
    expect_equal(result$iterations_completed, 0L)
    expect_equal(result$design_control, controls)
    expect_equal(result$model$n, 5L)
  }
})

test_that("invalid new design controls fail before any objective evaluations", {
  calls <- 0L
  objective <- function(...) {
    calls <<- calls + 1L
    0
  }
  for (control in list(
    list(sequence_method = "unknown"),
    list(sequence_method = "sfta", sfta_control = list(nrounds = 41L)),
    list(sequence_method = "random", sfta_control = list(nstarts = 2L))
  )) {
    expect_error(magp_bayes_optimize_from_scratch(
      FUN = objective, n_initial = 4, q = 3, n_iter = 0,
      design_control = c(list(sequence_maxit = 40, quantity_maxit = 2,
                              alignment_maxit = 2), control),
      fit_control = list(maxeval = 2), seed = 622, verbose = FALSE
    ))
    expect_identical(calls, 0L)
  }
})
