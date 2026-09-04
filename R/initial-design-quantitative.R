# Helpers and public functions for the quantitative portion of an initial
# design. Column-wise swaps preserve the Latin hypercube structure throughout
# the search.

#' Validate a quantitative design matrix
#'
#' @noRd
.magp_validate_quantity_design <- function(quantity, n = NULL, q = NULL,
                                           name = "quantity",
                                           require_latin = FALSE) {
  if (!is.matrix(quantity) && !is.data.frame(quantity)) {
    stop(name, " must be a numeric matrix or data frame", call. = FALSE)
  }
  quantity <- as.matrix(quantity)
  if (!is.numeric(quantity) || any(!is.finite(quantity))) {
    stop(name, " must contain only finite numeric values", call. = FALSE)
  }
  if (nrow(quantity) < 2L || ncol(quantity) < 1L) {
    stop(name, " must have at least two rows and one column", call. = FALSE)
  }
  if (!is.null(n) && nrow(quantity) != n) {
    stop(name, " must have ", n, " rows", call. = FALSE)
  }
  if (!is.null(q) && ncol(quantity) != q) {
    stop(name, " must have ", q, " columns", call. = FALSE)
  }
  if (any(quantity < 0 | quantity > 1)) {
    stop(name, " must lie in the interval [0, 1]", call. = FALSE)
  }

  if (require_latin) {
    n_runs <- nrow(quantity)
    upper <- 1 - .Machine$double.eps
    strata <- floor(pmin(quantity, upper) * n_runs)
    valid <- apply(strata, 2L, function(column) {
      identical(sort(as.integer(column)), 0:(n_runs - 1L))
    })
    if (any(!valid)) {
      bad <- which(!valid)
      shown <- paste(bad[seq_len(min(length(bad), 5L))], collapse = ", ")
      suffix <- if (length(bad) > 5L) ", ..." else ""
      stop(
        name, " must contain every Latin-hypercube stratum once in column(s): ",
        shown, suffix,
        call. = FALSE
      )
    }
  }

  storage.mode(quantity) <- "double"
  quantity
}

#' Calculate a stable inverse-distance p-norm
#'
#' @noRd
.magp_inverse_distance_norm <- function(distances, p) {
  if (!length(distances)) {
    stop("at least one pairwise distance is required", call. = FALSE)
  }
  if (any(!is.finite(distances)) || any(distances < 0)) {
    stop("pairwise distances must be finite and nonnegative", call. = FALSE)
  }
  if (any(distances == 0)) {
    return(Inf)
  }
  smallest <- min(distances)
  sum((smallest / distances)^p)^(1 / p) / smallest
}

#' Evaluate a quantitative Latin hypercube
#'
#' Measures the separation of design rows under Euclidean distance. The
#' criterion is an inverse-distance p-norm, so smaller values indicate a more
#' space-filling design.
#'
#' @param quantity Numeric matrix or data frame with values in `[0, 1]`.
#' @param p Positive whole-number exponent controlling emphasis on the
#'   shortest pairwise distances.
#'
#' @return One numeric criterion value. Smaller values are preferred. A design
#'   with duplicate rows has value `Inf`.
#'
#' @examples
#' quantity <- cbind(
#'   c(0.125, 0.375, 0.625, 0.875),
#'   c(0.625, 0.125, 0.875, 0.375)
#' )
#' magp_quantitative_criterion(quantity)
#'
#' @export
magp_quantitative_criterion <- function(quantity, p = 15L) {
  quantity <- .magp_validate_quantity_design(quantity)
  .magp2d_validate_scalar(p, "p", lower = 1, integer = TRUE)
  distances <- as.numeric(stats::dist(quantity, method = "euclidean"))
  .magp_inverse_distance_norm(distances, as.integer(p))
}

#' Generate a random Latin hypercube
#'
#' @noRd
.magp_random_latin_hypercube <- function(n, q) {
  strata <- vapply(
    seq_len(q),
    function(column) sample.int(n),
    integer(n)
  )
  quantity <- (strata - 0.5) / n
  colnames(quantity) <- paste0("quantity_", seq_len(q))
  quantity
}

#' Construct the quantitative portion of an initial design
#'
#' Uses simulated annealing to improve a Latin hypercube under a maximin-style
#' Euclidean-distance criterion. A candidate is created by swapping two values
#' within one column, which preserves the Latin hypercube property.
#'
#' @param n Number of design runs. Must be at least two.
#' @param q Number of quantitative variables. Must be positive.
#' @param p Positive whole-number exponent controlling emphasis on the
#'   shortest pairwise distances.
#' @param maxit Positive whole number of simulated-annealing iterations.
#' @param temp Positive initial temperature passed to [stats::optim()].
#' @param tmax Positive whole number of evaluations at each temperature.
#' @param initial Optional `n` by `q` Latin hypercube with values in `[0, 1]`.
#' @param seed Optional nonnegative whole-number seed. When supplied, the
#'   function restores the caller's random-number state before returning.
#'
#' @return An object of class `magp_quantitative_design`, containing the
#'   optimized Latin hypercube, the starting design, criterion values, minimum
#'   distances, and search settings.
#'
#' @references
#' Kirkpatrick, S., Gelatt, C. D., and Vecchi, M. P. (1983).
#' Optimization by Simulated Annealing. Science, 220, 671-680.
#' \doi{10.1126/science.220.4598.671}.
#'
#' @examples
#' design <- magp_quantitative_design(
#'   n = 8,
#'   q = 4,
#'   maxit = 500,
#'   seed = 1
#' )
#' design$quantity
#' design$criterion
#'
#' @export
magp_quantitative_design <- function(n, q, p = 15L, maxit = 10000L,
                                     temp = 0.01, tmax = 10L,
                                     initial = NULL, seed = NULL) {
  call <- match.call()
  .magp2d_validate_scalar(n, "n", lower = 2, integer = TRUE)
  .magp2d_validate_scalar(q, "q", lower = 1, integer = TRUE)
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

  output <- .magp_with_local_seed(seed, {
    if (is.null(initial)) {
      initial_quantity <- .magp_random_latin_hypercube(n, q)
    } else {
      initial_quantity <- .magp_validate_quantity_design(
        initial, n = n, q = q, name = "initial", require_latin = TRUE
      )
      colnames(initial_quantity) <- paste0("quantity_", seq_len(q))
    }
    initial_parameter <- as.vector(t(initial_quantity))

    objective <- function(parameter) {
      candidate <- matrix(parameter, nrow = n, ncol = q, byrow = TRUE)
      distances <- as.numeric(stats::dist(candidate, method = "euclidean"))
      .magp_inverse_distance_norm(distances, p)
    }
    neighbor <- function(parameter) {
      candidate <- matrix(parameter, nrow = n, ncol = q, byrow = TRUE)
      selected_column <- sample.int(q, 1L)
      selected_rows <- sample.int(n, 2L, replace = FALSE)
      temporary <- candidate[selected_rows[1L], selected_column]
      candidate[selected_rows[1L], selected_column] <-
        candidate[selected_rows[2L], selected_column]
      candidate[selected_rows[2L], selected_column] <- temporary
      as.vector(t(candidate))
    }

    fit <- stats::optim(
      par = initial_parameter,
      fn = objective,
      gr = neighbor,
      method = "SANN",
      control = list(maxit = maxit, temp = temp, tmax = tmax)
    )
    final_quantity <- matrix(fit$par, nrow = n, ncol = q, byrow = TRUE)
    initial_criterion <- objective(initial_parameter)
    final_criterion <- objective(as.vector(t(final_quantity)))
    if (final_criterion > initial_criterion) {
      final_quantity <- initial_quantity
      final_criterion <- initial_criterion
    }
    colnames(final_quantity) <- paste0("quantity_", seq_len(q))

    list(
      quantity = final_quantity,
      initial_quantity = initial_quantity,
      criterion = as.numeric(final_criterion),
      initial_criterion = as.numeric(initial_criterion),
      improvement = as.numeric(initial_criterion - final_criterion),
      minimum_distance = min(stats::dist(final_quantity)),
      initial_minimum_distance = min(stats::dist(initial_quantity)),
      n = n,
      q = q,
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
  class(output) <- "magp_quantitative_design"
  output
}

#' Print a quantitative initial design
#'
#' @param x A `magp_quantitative_design` object.
#' @param ... Additional arguments, currently unused.
#'
#' @return `x`, invisibly.
#'
#' @export
print.magp_quantitative_design <- function(x, ...) {
  cat("Quantitative initial design\n")
  cat("  Runs:", x$n, "\n")
  cat("  Variables:", x$q, "\n")
  cat("  Method:", x$method, "\n")
  cat("  Initial criterion:", format(x$initial_criterion, digits = 7), "\n")
  cat("  Final criterion:", format(x$criterion, digits = 7), "\n")
  cat("  Minimum distance:", format(x$minimum_distance, digits = 7), "\n")
  invisible(x)
}
