reference_sequence_criterion <- function(sequence, pair_weight = 0.2,
                                         space_weight = 0.8, p = 15L) {
  n <- nrow(sequence)
  q <- ncol(sequence)
  orders <- t(vapply(seq_len(n), function(i) {
    order(sequence[i, ])
  }, integer(q)))

  pair_counts <- matrix(0, nrow = q, ncol = q)
  for (i in seq_len(n)) {
    for (j in seq_len(q - 1L)) {
      pair_counts[orders[i, j], orders[i, j + 1L]] <-
        pair_counts[orders[i, j], orders[i, j + 1L]] + 1
    }
  }
  pair_term <- sum((pair_counts[row(pair_counts) != col(pair_counts)] + 1)^(-p))

  row_pairs <- utils::combn(n, 2L)
  hamming <- apply(row_pairs, 2L, function(index) {
    sum(orders[index[1L], ] != orders[index[2L], ])
  })
  space_term <- sum((hamming + 1)^(-p))

  total_weight <- pair_weight + space_weight
  (
    pair_weight / total_weight * pair_term +
      space_weight / total_weight * space_term
  )^(1 / p)
}

test_that("sequence criterion matches an independent calculation", {
  sequence <- rbind(
    c(1, 2, 3, 4),
    c(3, 1, 4, 2),
    c(2, 4, 1, 3),
    c(4, 3, 2, 1)
  )

  expect_equal(
    magp_sequence_criterion(sequence),
    reference_sequence_criterion(sequence),
    tolerance = 1e-14
  )
})

test_that("balanced and separated sequence designs score better", {
  balanced <- rbind(
    c(1, 2, 3, 4),
    c(3, 1, 4, 2),
    c(2, 4, 1, 3),
    c(4, 3, 2, 1)
  )
  repeated <- matrix(rep(1:4, 4L), nrow = 4L, byrow = TRUE)

  expect_lt(
    magp_sequence_criterion(balanced),
    magp_sequence_criterion(repeated)
  )
})

test_that("simulated annealing returns valid and reproducible sequences", {
  first <- magp_sequence_design(
    n = 8,
    q = 4,
    maxit = 500,
    temp = 0.1,
    seed = 317
  )
  second <- magp_sequence_design(
    n = 8,
    q = 4,
    maxit = 500,
    temp = 0.1,
    seed = 317
  )

  expect_s3_class(first, "magp_sequence_design")
  expect_equal(first$sequence, second$sequence)
  expect_equal(first$criterion, second$criterion)
  expect_equal(dim(first$sequence), c(8L, 4L))
  expect_true(all(apply(first$sequence, 1L, function(row) {
    identical(sort(as.integer(row)), 1:4)
  })))
  expect_lte(first$criterion, first$initial_criterion)
  expect_equal(
    first$criterion,
    magp_sequence_criterion(first$sequence),
    tolerance = 1e-14
  )
})

test_that("sequence search works with different component counts", {
  for (q in c(3L, 5L, 6L)) {
    design <- magp_sequence_design(
      n = 2L * q,
      q = q,
      maxit = 100,
      temp = 0.1,
      seed = 100 + q
    )

    expect_equal(dim(design$sequence), c(2L * q, q))
    expect_true(all(apply(design$sequence, 1L, function(row) {
      identical(sort(as.integer(row)), seq_len(q))
    })))
    expect_lte(design$criterion, design$initial_criterion)
  }
})

test_that("a supplied seed does not change the caller's random state", {
  set.seed(99)
  state_before <- .Random.seed
  invisible(magp_sequence_design(n = 5, q = 4, maxit = 50, seed = 12))
  expect_identical(.Random.seed, state_before)
})

test_that("sequence design inputs are checked before optimization", {
  expect_error(magp_sequence_design(n = 1, q = 4), "n")
  expect_error(magp_sequence_design(n = 4, q = 2), "q")
  expect_error(magp_sequence_design(n = 4, q = 4, p = 0), "p")
  expect_error(magp_sequence_design(n = 4, q = 4, maxit = 0), "maxit")
  expect_error(magp_sequence_design(n = 4, q = 4, temp = 0), "temp")
  expect_error(
    magp_sequence_design(
      n = 4, q = 4, pair_weight = 0, space_weight = 0
    ),
    "at least one"
  )

  invalid <- matrix(rep(1:4, 4L), nrow = 4L, byrow = TRUE)
  invalid[1, ] <- c(1, 1, 3, 4)
  expect_error(
    magp_sequence_design(n = 4, q = 4, initial = invalid),
    "permutation"
  )
  expect_error(
    magp_sequence_criterion(matrix(1:8, nrow = 2L)),
    "permutation"
  )
})

test_that("print method reports the main search results", {
  design <- magp_sequence_design(n = 4, q = 4, maxit = 20, seed = 5)
  expect_output(print(design), "Sequence initial design")
  expect_output(print(design), "Final criterion")
})
