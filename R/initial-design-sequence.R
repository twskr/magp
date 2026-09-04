# Helpers and public functions for the sequence portion of a quantitative-
# sequence initial design. The quantitative portion and joint alignment are
# intentionally handled separately so each stage can be assessed on its own.

#' Validate a matrix of sequence positions
#'
#' @noRd
.magp_validate_sequence_positions <- function(sequence, n = NULL, q = NULL,
                                               name = "sequence") {
  if (!is.matrix(sequence) && !is.data.frame(sequence)) {
    stop(name, " must be a numeric matrix or data frame", call. = FALSE)
  }
  sequence <- as.matrix(sequence)
  if (!is.numeric(sequence) || any(!is.finite(sequence))) {
    stop(name, " must contain only finite numeric values", call. = FALSE)
  }
  if (!is.null(n) && nrow(sequence) != n) {
    stop(name, " must have ", n, " rows", call. = FALSE)
  }
  if (!is.null(q) && ncol(sequence) != q) {
    stop(name, " must have ", q, " columns", call. = FALSE)
  }

  q <- ncol(sequence)
  if (q < 3L) {
    stop(name, " must contain at least three sequence columns", call. = FALSE)
  }
  valid <- apply(sequence, 1L, function(row) {
    all(row == floor(row)) &&
      identical(sort(as.integer(row)), seq_len(q))
  })
  if (any(!valid)) {
    bad <- which(!valid)
    shown <- paste(bad[seq_len(min(length(bad), 5L))], collapse = ", ")
    suffix <- if (length(bad) > 5L) ", ..." else ""
    stop(
      name, " must be a permutation of 1:q in every row; invalid row(s): ",
      shown, suffix,
      call. = FALSE
    )
  }
  storage.mode(sequence) <- "integer"
  sequence
}

#' Normalize the two sequence-design weights
#'
#' @noRd
.magp_sequence_weights <- function(pair_weight, space_weight) {
  .magp2d_validate_scalar(pair_weight, "pair_weight", lower = 0)
  .magp2d_validate_scalar(space_weight, "space_weight", lower = 0)
  total <- pair_weight + space_weight
  if (total <= 0) {
    stop(
      "at least one of pair_weight and space_weight must be positive",
      call. = FALSE
    )
  }
  c(pair_balance = pair_weight / total, space_filling = space_weight / total)
}

#' Convert component positions to components in execution order
#'
#' @noRd
.magp_positions_to_orders <- function(sequence) {
  q <- ncol(sequence)
  t(vapply(
    seq_len(nrow(sequence)),
    function(i) order(sequence[i, ]),
    integer(q)
  ))
}

#' Convert components in execution order to component positions
#'
#' @noRd
.magp_orders_to_positions <- function(orders) {
  q <- ncol(orders)
  positions <- t(vapply(
    seq_len(nrow(orders)),
    function(i) order(orders[i, ]),
    integer(q)
  ))
  colnames(positions) <- paste0("sequence_", seq_len(q))
  positions
}

#' Evaluate the sequence-design criterion in execution-order form
#'
#' @noRd
.magp_sequence_criterion_orders <- function(orders, weights, p) {
  n <- nrow(orders)
  q <- ncol(orders)

  from <- as.vector(orders[, -q, drop = FALSE])
  to <- as.vector(orders[, -1L, drop = FALSE])
  pair_index <- (from - 1L) * q + to
  pair_counts <- matrix(
    tabulate(pair_index, nbins = q * q),
    nrow = q,
    byrow = TRUE
  )
  pair_term <- sum((pair_counts[row(pair_counts) != col(pair_counts)] + 1)^(-p))

  row_pairs <- utils::combn(n, 2L)
  hamming <- vapply(seq_len(ncol(row_pairs)), function(i) {
    sum(orders[row_pairs[1L, i], ] != orders[row_pairs[2L, i], ])
  }, numeric(1))
  space_term <- sum((hamming + 1)^(-p))

  (
    weights[["pair_balance"]] * pair_term +
      weights[["space_filling"]] * space_term
  )^(1 / p)
}

#' Evaluate a sequence initial design
#'
#' Measures two properties of a sequence design: how evenly ordered adjacent
#' component pairs are represented and how well separated the design rows are
#' under Hamming distance. Smaller values indicate a better design.
#'
#' Each row uses the same sequence format as the fitting functions. The value
#' in column `j` is the position assigned to component `j`, and every row must
#' be a permutation of `1:q`.
#'
#' @param sequence Numeric matrix or data frame. Every row must be a
#'   permutation of `1:q`.
#' @param pair_weight Nonnegative weight for ordered adjacent-pair balance.
#' @param space_weight Nonnegative weight for Hamming-distance space filling.
#' @param p Positive whole-number exponent controlling emphasis on the weakest
#'   pair counts and the smallest distances.
#'
#' @return One numeric criterion value. Smaller values are preferred.
#'
#' @references
#' Xiao, Q., Wang, Y., Mandal, A., and Deng, X. (2024).
#' Modeling and Active Learning for Experiments with Quantitative-Sequence
#' Factors. Journal of the American Statistical Association.
#' \doi{10.1080/01621459.2022.2123335}.
#'
#' @examples
#' sequence <- rbind(
#'   c(1, 2, 3, 4),
#'   c(3, 1, 4, 2),
#'   c(2, 4, 1, 3),
#'   c(4, 3, 2, 1)
#' )
#' magp_sequence_criterion(sequence)
#'
#' @export
magp_sequence_criterion <- function(sequence, pair_weight = 0.2,
                                    space_weight = 0.8, p = 15L) {
  sequence <- .magp_validate_sequence_positions(sequence)
  .magp2d_validate_scalar(p, "p", lower = 1, integer = TRUE)
  weights <- .magp_sequence_weights(pair_weight, space_weight)
  orders <- .magp_positions_to_orders(sequence)
  .magp_sequence_criterion_orders(orders, weights, as.integer(p))
}

#' Generate a random sequence design
#'
#' @noRd
.magp_random_sequence_positions <- function(n, q) {
  sequence <- matrix(NA_integer_, nrow = n, ncol = q)
  keys <- new.env(hash = TRUE, parent = emptyenv())
  unique_rows_possible <- lgamma(q + 1) >= log(n)

  for (i in seq_len(n)) {
    candidate <- sample.int(q)
    if (unique_rows_possible) {
      attempt <- 0L
      key <- paste(candidate, collapse = ",")
      while (exists(key, envir = keys, inherits = FALSE) && attempt < 10000L) {
        candidate <- sample.int(q)
        key <- paste(candidate, collapse = ",")
        attempt <- attempt + 1L
      }
      if (exists(key, envir = keys, inherits = FALSE)) {
        stop("could not generate a unique starting sequence design", call. = FALSE)
      }
      assign(key, TRUE, envir = keys)
    }
    sequence[i, ] <- candidate
  }
  colnames(sequence) <- paste0("sequence_", seq_len(q))
  sequence
}

#' Evaluate code with a temporary random-number seed
#'
#' @noRd
.magp_with_local_seed <- function(seed, code) {
  if (is.null(seed)) {
    return(force(code))
  }
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  force(code)
}

#' Construct the sequence portion of an initial design
#'
#' Uses simulated annealing to search for a sequence design with balanced
#' ordered adjacent pairs and well-separated rows. A neighbor is generated by
#' selecting one design row and exchanging two component positions, so every
#' candidate remains a valid permutation.
#'
#' This function constructs only the sequence portion of a quantitative-
#' sequence initial design. Use [magp_initial_design()] to combine it with a
#' quantitative Latin hypercube and improve the pairing of the two portions.
#'
#' @param n Number of design runs. Must be at least two.
#' @param q Number of components. Must be at least three.
#' @param pair_weight Nonnegative weight for ordered adjacent-pair balance.
#' @param space_weight Nonnegative weight for Hamming-distance space filling.
#' @param p Positive whole-number exponent controlling emphasis on the weakest
#'   pair counts and the smallest distances.
#' @param maxit Positive whole number of simulated-annealing iterations.
#' @param temp Positive initial temperature passed to [stats::optim()]. The
#'   default is scaled to the sequence-design criterion used here.
#' @param tmax Positive whole number of evaluations at each temperature.
#' @param initial Optional `n` by `q` matrix of component positions. Every row
#'   must be a permutation of `1:q`.
#' @param seed Optional nonnegative whole-number seed. When supplied, the
#'   function restores the caller's random-number state before returning.
#'
#' @return An object of class `magp_sequence_design`, containing the optimized
#'   sequence matrix, the starting matrix, criterion values, and search
#'   settings. The sequence matrix is ready to use as the sequence half of a
#'   `magp` input matrix.
#'
#' @references
#' Kirkpatrick, S., Gelatt, C. D., and Vecchi, M. P. (1983).
#' Optimization by Simulated Annealing. Science, 220, 671-680.
#' \doi{10.1126/science.220.4598.671}.
#'
#' Xiao, Q., Wang, Y., Mandal, A., and Deng, X. (2024).
#' Modeling and Active Learning for Experiments with Quantitative-Sequence
#' Factors. Journal of the American Statistical Association.
#' \doi{10.1080/01621459.2022.2123335}.
#'
#' @examples
#' design <- magp_sequence_design(
#'   n = 8,
#'   q = 4,
#'   maxit = 500,
#'   seed = 1
#' )
#' design$sequence
#' design$criterion
#'
#' @export
magp_sequence_design <- function(n, q, pair_weight = 0.2,
                                 space_weight = 0.8, p = 15L,
                                 maxit = 10000L, temp = 0.1, tmax = 10L,
                                 initial = NULL, seed = NULL) {
  call <- match.call()
  .magp2d_validate_scalar(n, "n", lower = 2, integer = TRUE)
  .magp2d_validate_scalar(q, "q", lower = 3, integer = TRUE)
  .magp2d_validate_scalar(p, "p", lower = 1, integer = TRUE)
  .magp2d_validate_scalar(maxit, "maxit", lower = 1, integer = TRUE)
  .magp2d_validate_scalar(temp, "temp", lower = 0, strict_lower = TRUE)
  .magp2d_validate_scalar(tmax, "tmax", lower = 1, integer = TRUE)
  .magp2d_validate_scalar(
    seed, "seed", lower = 0, upper = .Machine$integer.max,
    integer = TRUE, allow_null = TRUE
  )

  n <- as.integer(n)
  q <- as.integer(q)
  p <- as.integer(p)
  maxit <- as.integer(maxit)
  tmax <- as.integer(tmax)
  weights <- .magp_sequence_weights(pair_weight, space_weight)

  output <- .magp_with_local_seed(seed, {
    if (is.null(initial)) {
      initial_sequence <- .magp_random_sequence_positions(n, q)
    } else {
      initial_sequence <- .magp_validate_sequence_positions(
        initial, n = n, q = q, name = "initial"
      )
      colnames(initial_sequence) <- paste0("sequence_", seq_len(q))
    }
    initial_orders <- .magp_positions_to_orders(initial_sequence)
    initial_parameter <- as.vector(t(initial_orders))

    objective <- function(parameter) {
      orders <- matrix(
        as.integer(round(parameter)), nrow = n, ncol = q, byrow = TRUE
      )
      .magp_sequence_criterion_orders(orders, weights, p)
    }
    neighbor <- function(parameter) {
      orders <- matrix(
        as.integer(round(parameter)), nrow = n, ncol = q, byrow = TRUE
      )
      selected_row <- sample.int(n, 1L)
      selected_columns <- sample.int(q, 2L, replace = FALSE)
      temporary <- orders[selected_row, selected_columns[1L]]
      orders[selected_row, selected_columns[1L]] <-
        orders[selected_row, selected_columns[2L]]
      orders[selected_row, selected_columns[2L]] <- temporary
      as.vector(t(orders))
    }

    fit <- stats::optim(
      par = initial_parameter,
      fn = objective,
      gr = neighbor,
      method = "SANN",
      control = list(maxit = maxit, temp = temp, tmax = tmax)
    )
    final_orders <- matrix(
      as.integer(round(fit$par)), nrow = n, ncol = q, byrow = TRUE
    )
    initial_criterion <- objective(initial_parameter)
    final_criterion <- objective(as.vector(t(final_orders)))
    if (final_criterion > initial_criterion) {
      final_orders <- initial_orders
      final_criterion <- initial_criterion
    }

    list(
      sequence = .magp_orders_to_positions(final_orders),
      initial_sequence = initial_sequence,
      criterion = as.numeric(final_criterion),
      initial_criterion = as.numeric(initial_criterion),
      improvement = as.numeric(initial_criterion - final_criterion),
      n = n,
      q = q,
      weights = weights,
      p = p,
      method = "simulated annealing",
      maxit = maxit,
      temp = temp,
      tmax = tmax,
      seed = seed,
      evaluations = unname(fit$counts),
      convergence = fit$convergence
    )
  })
  output$call <- call
  class(output) <- "magp_sequence_design"
  output
}

#' Print a sequence initial design
#'
#' @param x A fitted `magp_sequence_design` object.
#' @param ... Additional arguments, currently unused.
#'
#' @return `x`, invisibly.
#'
#' @export
print.magp_sequence_design <- function(x, ...) {
  cat("Sequence initial design\n")
  cat("  Runs:", x$n, "\n")
  cat("  Components:", x$q, "\n")
  cat("  Method:", x$method, "\n")
  cat("  Initial criterion:", format(x$initial_criterion, digits = 7), "\n")
  cat("  Final criterion:", format(x$criterion, digits = 7), "\n")
  invisible(x)
}
