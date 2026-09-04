latin_hypercube_is_valid <- function(quantity) {
  n <- nrow(quantity)
  strata <- floor(pmin(quantity, 1 - .Machine$double.eps) * n)
  all(apply(strata, 2L, function(column) {
    identical(sort(as.integer(column)), 0:(n - 1L))
  }))
}

reference_quantitative_criterion <- function(quantity, p = 15L) {
  distances <- as.numeric(stats::dist(quantity))
  sum(distances^(-p))^(1 / p)
}

reference_joint_criterion <- function(quantity, sequence,
                                      quantity_weight = 0.5,
                                      sequence_weight = 0.5, p = 15L) {
  total_weight <- quantity_weight + sequence_weight
  quantity_weight <- quantity_weight / total_weight
  sequence_weight <- sequence_weight / total_weight
  pairs <- utils::combn(nrow(quantity), 2L)
  combined <- apply(pairs, 2L, function(index) {
    quantity_distance <- sqrt(sum(
      (quantity[index[1L], ] - quantity[index[2L], ])^2
    ))
    sequence_distance <- sum(
      sequence[index[1L], ] != sequence[index[2L], ]
    )
    quantity_weight * quantity_distance +
      sequence_weight * sequence_distance + 1
  })
  sum(combined^(-p))^(1 / p)
}

test_that("quantitative criterion matches an independent calculation", {
  quantity <- cbind(
    c(0.125, 0.375, 0.625, 0.875),
    c(0.625, 0.125, 0.875, 0.375)
  )

  expect_equal(
    magp_quantitative_criterion(quantity),
    reference_quantitative_criterion(quantity),
    tolerance = 1e-14
  )
  expect_identical(
    magp_quantitative_criterion(rbind(c(0, 0), c(0, 0))),
    Inf
  )
})

test_that("quantitative search preserves the Latin hypercube", {
  first <- magp_quantitative_design(
    n = 8, q = 4, maxit = 300, temp = 0.1, seed = 27
  )
  second <- magp_quantitative_design(
    n = 8, q = 4, maxit = 300, temp = 0.1, seed = 27
  )

  expect_s3_class(first, "magp_quantitative_design")
  expect_equal(first$quantity, second$quantity)
  expect_equal(first$criterion, second$criterion)
  expect_equal(dim(first$quantity), c(8L, 4L))
  expect_true(latin_hypercube_is_valid(first$quantity))
  expect_lte(first$criterion, first$initial_criterion)
  expect_equal(
    first$criterion,
    magp_quantitative_criterion(first$quantity),
    tolerance = 1e-14
  )
})

test_that("quantitative search works with different dimensions", {
  for (q in c(1L, 3L, 5L)) {
    design <- magp_quantitative_design(
      n = 2L * q + 2L,
      q = q,
      maxit = 100,
      temp = 0.1,
      seed = 200 + q
    )
    expect_equal(dim(design$quantity), c(2L * q + 2L, q))
    expect_true(latin_hypercube_is_valid(design$quantity))
    expect_lte(design$criterion, design$initial_criterion)
  }
})

test_that("quantitative inputs are checked", {
  expect_error(magp_quantitative_design(n = 1, q = 3), "n")
  expect_error(magp_quantitative_design(n = 4, q = 0), "q")
  expect_error(magp_quantitative_design(n = 4, q = 3, p = 0), "p")
  expect_error(magp_quantitative_design(n = 4, q = 3, maxit = 0), "maxit")
  expect_error(magp_quantitative_design(n = 4, q = 3, temp = 0), "temp")
  expect_error(
    magp_quantitative_criterion(matrix(c(0, 0.5, 1.2, 0.4), ncol = 2)),
    "\\[0, 1\\]"
  )

  invalid_latin <- matrix(
    c(0.125, 0.125, 0.625, 0.875, 0.125, 0.375, 0.625, 0.875),
    nrow = 4L
  )
  expect_error(
    magp_quantitative_design(
      n = 4, q = 2, initial = invalid_latin, maxit = 10
    ),
    "Latin-hypercube"
  )
})

test_that("joint criterion matches an independent calculation", {
  quantity <- cbind(
    c(0.125, 0.375, 0.625, 0.875),
    c(0.625, 0.125, 0.875, 0.375),
    c(0.375, 0.875, 0.125, 0.625),
    c(0.875, 0.625, 0.375, 0.125)
  )
  sequence <- rbind(
    c(1, 2, 3, 4),
    c(3, 1, 4, 2),
    c(2, 4, 1, 3),
    c(4, 3, 2, 1)
  )

  expect_equal(
    magp_joint_criterion(quantity, sequence),
    reference_joint_criterion(quantity, sequence),
    tolerance = 1e-14
  )
  expect_equal(
    magp_joint_criterion(
      quantity, sequence, quantity_weight = 2, sequence_weight = 2
    ),
    magp_joint_criterion(quantity, sequence),
    tolerance = 1e-14
  )
  expect_error(
    magp_joint_criterion(
      quantity, sequence, quantity_weight = 0, sequence_weight = 0
    ),
    "at least one"
  )
})

test_that("complete initial design preserves both component designs", {
  first <- magp_initial_design(
    n = 8,
    q = 4,
    sequence_maxit = 300,
    quantity_maxit = 300,
    alignment_maxit = 300,
    seed = 42
  )
  second <- magp_initial_design(
    n = 8,
    q = 4,
    sequence_maxit = 300,
    quantity_maxit = 300,
    alignment_maxit = 300,
    seed = 42
  )

  expect_s3_class(first, "magp_initial_design")
  expect_equal(first$design, second$design)
  expect_equal(dim(first$design), c(8L, 8L))
  expect_named(
    as.data.frame(first$design),
    c(paste0("quantity_", 1:4), paste0("sequence_", 1:4))
  )
  expect_true(latin_hypercube_is_valid(first$quantity))
  expect_true(all(apply(first$sequence, 1L, function(row) {
    identical(sort(as.integer(row)), 1:4)
  })))
  expect_equal(first$sequence, first$sequence_search$sequence)
  expect_equal(
    magp_quantitative_criterion(first$quantity),
    first$quantity_search$criterion,
    tolerance = 1e-14
  )
  expect_lte(
    first$criteria[["joint"]],
    first$criteria[["joint_before_alignment"]]
  )
  expect_equal(
    first$criteria[["joint"]],
    magp_joint_criterion(first$quantity, first$sequence),
    tolerance = 1e-14
  )

  prepared <- magp:::.magp2d_prepare_fit_inputs(
    first$design, y = seq_len(8L), q = 4L
  )
  expect_equal(dim(prepared$X), c(8L, 8L))
  expect_false(any(prepared$quantity_scaling$applied))
})

test_that("complete design uses local seeds", {
  set.seed(991)
  state_before <- .Random.seed
  invisible(magp_initial_design(
    n = 6,
    q = 3,
    sequence_maxit = 50,
    quantity_maxit = 50,
    alignment_maxit = 50,
    seed = 8
  ))
  expect_identical(.Random.seed, state_before)
})

test_that("complete design inputs and print method are clear", {
  expect_error(magp_initial_design(n = 4, q = 2), "q")
  expect_error(
    magp_initial_design(
      n = 4, q = 3, alignment_maxit = 0,
      sequence_maxit = 10, quantity_maxit = 10
    ),
    "alignment_maxit"
  )
  expect_error(
    magp_joint_criterion(
      matrix(runif(12), nrow = 4),
      matrix(rep(1:4, 4), nrow = 4, byrow = TRUE)
    ),
    "columns"
  )

  design <- magp_initial_design(
    n = 4,
    q = 3,
    sequence_maxit = 20,
    quantity_maxit = 20,
    alignment_maxit = 20,
    seed = 3
  )
  expect_output(print(design), "Quantitative-sequence initial design")
  expect_output(print(design$quantity_search), "Minimum distance")
})
