# Functions that combine the quantitative and sequence portions of an initial
# design. The row-alignment search preserves both component designs.

#' Normalize the weights used in the joint design criterion
#'
#' @noRd
.magp_joint_weights <- function(quantity_weight, sequence_weight) {
  .magp2d_validate_scalar(quantity_weight, "quantity_weight", lower = 0)
  .magp2d_validate_scalar(sequence_weight, "sequence_weight", lower = 0)
  total <- quantity_weight + sequence_weight
  if (total <= 0) {
    stop(
      "at least one of quantity_weight and sequence_weight must be positive",
      call. = FALSE
    )
  }
  c(
    quantity = quantity_weight / total,
    sequence = sequence_weight / total
  )
}

#' Evaluate the joint design criterion after input validation
#'
#' @noRd
.magp_joint_criterion_values <- function(quantity, sequence, weights, p) {
  n <- nrow(quantity)
  row_pairs <- utils::combn(n, 2L)
  quantity_distance <- vapply(seq_len(ncol(row_pairs)), function(i) {
    difference <-
      quantity[row_pairs[1L, i], ] - quantity[row_pairs[2L, i], ]
    sqrt(sum(difference^2))
  }, numeric(1))
  sequence_distance <- vapply(seq_len(ncol(row_pairs)), function(i) {
    sum(sequence[row_pairs[1L, i], ] != sequence[row_pairs[2L, i], ])
  }, numeric(1))
  combined_distance <-
    weights[["quantity"]] * quantity_distance +
      weights[["sequence"]] * sequence_distance + 1
  .magp_inverse_distance_norm(combined_distance, p)
}

#' Evaluate a complete quantitative-sequence initial design
#'
#' Combines Euclidean distance for the quantitative portion with Hamming
#' distance for the sequence portion. Smaller values indicate better joint
#' separation of the design runs.
#'
#' @param quantity Numeric matrix or data frame with values in `[0, 1]`.
#' @param sequence Numeric matrix or data frame with the same dimensions as
#'   `quantity`. Every row must be a permutation of `1:q`.
#' @param quantity_weight Nonnegative weight for quantitative distance.
#' @param sequence_weight Nonnegative weight for sequence distance.
#' @param p Positive whole-number exponent controlling emphasis on the least
#'   separated pairs of runs.
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
#' quantity <- cbind(
#'   c(0.125, 0.375, 0.625, 0.875),
#'   c(0.625, 0.125, 0.875, 0.375),
#'   c(0.375, 0.875, 0.125, 0.625),
#'   c(0.875, 0.625, 0.375, 0.125)
#' )
#' sequence <- rbind(
#'   c(1, 2, 3, 4),
#'   c(3, 1, 4, 2),
#'   c(2, 4, 1, 3),
#'   c(4, 3, 2, 1)
#' )
#' magp_joint_criterion(quantity, sequence)
#'
#' @export
magp_joint_criterion <- function(quantity, sequence, quantity_weight = 0.5,
                                 sequence_weight = 0.5, p = 15L) {
  quantity <- .magp_validate_quantity_design(quantity)
  sequence <- .magp_validate_sequence_positions(
    sequence, n = nrow(quantity), q = ncol(quantity)
  )
  .magp2d_validate_scalar(p, "p", lower = 1, integer = TRUE)
  weights <- .magp_joint_weights(quantity_weight, sequence_weight)
  .magp_joint_criterion_values(quantity, sequence, weights, as.integer(p))
}

#' Offset a seed without exceeding R's integer range
#'
#' @noRd
.magp_design_seed <- function(seed, offset) {
  if (is.null(seed)) {
    return(NULL)
  }
  as.integer((as.double(seed) + offset) %% (.Machine$integer.max + 1))
}

#' Align fixed quantitative and sequence designs
#'
#' @noRd
.magp_align_initial_design <- function(quantity, sequence, weights, p,
                                       maxit, temp, tmax, seed) {
  n <- nrow(quantity)
  q <- ncol(quantity)
  row_pairs <- utils::combn(n, 2L)
  quantity_distances <- as.matrix(stats::dist(quantity, method = "euclidean"))
  sequence_distances <- vapply(seq_len(ncol(row_pairs)), function(i) {
    sum(sequence[row_pairs[1L, i], ] != sequence[row_pairs[2L, i], ])
  }, numeric(1))

  .magp_with_local_seed(seed, {
    initial_order <- seq_len(n)
    objective <- function(ordering) {
      ordering <- as.integer(round(ordering))
      quantity_distance <- quantity_distances[cbind(
        ordering[row_pairs[1L, ]], ordering[row_pairs[2L, ]]
      )]
      combined_distance <-
        weights[["quantity"]] * quantity_distance +
          weights[["sequence"]] * sequence_distances + 1
      .magp_inverse_distance_norm(combined_distance, p)
    }
    neighbor <- function(ordering) {
      ordering <- as.integer(round(ordering))
      selected <- sample.int(n, 2L, replace = FALSE)
      temporary <- ordering[selected[1L]]
      ordering[selected[1L]] <- ordering[selected[2L]]
      ordering[selected[2L]] <- temporary
      ordering
    }

    fit <- stats::optim(
      par = initial_order,
      fn = objective,
      gr = neighbor,
      method = "SANN",
      control = list(maxit = maxit, temp = temp, tmax = tmax)
    )
    final_order <- as.integer(round(fit$par))
    initial_criterion <- objective(initial_order)
    final_criterion <- objective(final_order)
    if (final_criterion > initial_criterion) {
      final_order <- initial_order
      final_criterion <- initial_criterion
    }

    list(
      quantity = quantity[final_order, , drop = FALSE],
      sequence = sequence,
      row_order = final_order,
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
}

#' Construct a quantitative-sequence initial design
#'
#' Builds an initial design in three steps. It first generates the sequence
#' permutations, then constructs a maximin-style Latin hypercube for the
#' quantitative levels, and finally aligns the two fixed designs by permuting
#' whole quantitative rows. The final alignment preserves both the Latin
#' hypercube and every sequence permutation.
#'
#' @param n Number of design runs. Must be at least two.
#' @param q Number of components. Must be at least three.
#' @param pair_weight Nonnegative weight for ordered adjacent-pair balance in
#'   the sequence search.
#' @param sequence_space_weight Nonnegative weight for Hamming-distance space
#'   filling in the sequence search.
#' @param quantity_weight Nonnegative quantitative-distance weight in the
#'   joint criterion.
#' @param sequence_weight Nonnegative sequence-distance weight in the joint
#'   criterion.
#' @param p Positive whole-number exponent used in all three criteria.
#' @param sequence_maxit Positive whole number of sequence-search iterations.
#' @param quantity_maxit Positive whole number of quantitative-search
#'   iterations.
#' @param alignment_maxit Positive whole number of row-alignment iterations.
#' @param sequence_temp Positive initial temperature for the sequence search.
#' @param quantity_temp Positive initial temperature for the quantitative
#'   search.
#' @param alignment_temp Positive initial temperature for the alignment search.
#' @param tmax Positive whole number of evaluations at each temperature.
#' @param initial_sequence Optional `n` by `q` matrix of sequence positions.
#' @param initial_quantity Optional `n` by `q` Latin hypercube in `[0, 1]`.
#' @param seed Optional nonnegative whole-number seed. Separate deterministic
#'   seeds are used for the three stages, and the caller's random-number state
#'   is preserved.
#' @param sequence_method Method for the sequence portion. Choose `"random"`,
#'   `"sfta"` for space-filling threshold accepting, or `"sann"` for simulated
#'   annealing. The default is `"sann"` for backward compatibility.
#' @param sfta_control Named list of SFTA settings, passed to
#'   [magp_sequence_design()]. Used only for `sequence_method = "sfta"`.
#'
#' @details The sequence method does not change how quantitative levels are
#'   generated. With a fixed `seed`, all three methods use the same quantitative
#'   design before the final row-alignment step. The alignment can change its
#'   row order but not its values or Latin-hypercube structure.
#'
#' @return An object of class `magp_initial_design`. Its `design` element is an
#'   `n` by `2*q` matrix ready to use as the input to a MAGP fitting function.
#'   The first `q` columns contain quantitative levels and the last `q` columns
#'   contain sequence positions. The component searches and their criterion
#'   values are retained for inspection.
#'
#' @references
#' Xiao, Q., Wang, Y., Mandal, A., and Deng, X. (2024).
#' Modeling and Active Learning for Experiments with Quantitative-Sequence
#' Factors. Journal of the American Statistical Association.
#' \doi{10.1080/01621459.2022.2123335}.
#'
#' @examples
#' design <- magp_initial_design(
#'   n = 8,
#'   q = 4,
#'   sequence_maxit = 500,
#'   quantity_maxit = 500,
#'   alignment_maxit = 500,
#'   seed = 1
#' )
#' design$design
#' design$criteria
#'
#' random <- magp_initial_design(
#'   8, 4, sequence_method = "random", seed = 1,
#'   quantity_maxit = 100, alignment_maxit = 100
#' )
#' sfta <- magp_initial_design(
#'   8, 4, sequence_method = "sfta", seed = 1, sequence_maxit = 200,
#'   quantity_maxit = 100, alignment_maxit = 100,
#'   sfta_control = list(nstarts = 2, ncalibrate = 50)
#' )
#' rbind(random = random$criteria, sfta = sfta$criteria)
#'
#' @export
magp_initial_design <- function(
    n, q, pair_weight = 0.2, sequence_space_weight = 0.8,
    quantity_weight = 0.5, sequence_weight = 0.5, p = 15L,
    sequence_maxit = 10000L, quantity_maxit = 10000L,
    alignment_maxit = 10000L, sequence_temp = 0.1,
    quantity_temp = 0.01, alignment_temp = 0.001, tmax = 10L,
    initial_sequence = NULL, initial_quantity = NULL, seed = NULL,
    sequence_method = c("sann", "random", "sfta"), sfta_control = list()) {
  call <- match.call()
  sequence_method <- match.arg(sequence_method)
  .magp2d_validate_scalar(n, "n", lower = 2, integer = TRUE)
  .magp2d_validate_scalar(q, "q", lower = 3, integer = TRUE)
  .magp2d_validate_scalar(p, "p", lower = 1, integer = TRUE)
  .magp2d_validate_scalar(
    alignment_maxit, "alignment_maxit", lower = 1, integer = TRUE
  )
  .magp2d_validate_scalar(
    alignment_temp, "alignment_temp", lower = 0, strict_lower = TRUE
  )
  .magp2d_validate_scalar(tmax, "tmax", lower = 1, integer = TRUE)
  .magp2d_validate_scalar(
    seed, "seed", lower = 0, upper = .Machine$integer.max,
    integer = TRUE, allow_null = TRUE
  )

  n <- as.integer(n)
  q <- as.integer(q)
  p <- as.integer(p)
  alignment_maxit <- as.integer(alignment_maxit)
  tmax <- as.integer(tmax)
  joint_weights <- .magp_joint_weights(quantity_weight, sequence_weight)

  sequence_search <- magp_sequence_design(
    n = n,
    q = q,
    pair_weight = pair_weight,
    space_weight = sequence_space_weight,
    p = p,
    maxit = sequence_maxit,
    temp = sequence_temp,
    tmax = tmax,
    initial = initial_sequence,
    seed = .magp_design_seed(seed, 0),
    method = sequence_method,
    sfta_control = sfta_control
  )
  quantity_search <- magp_quantitative_design(
    n = n,
    q = q,
    p = p,
    maxit = quantity_maxit,
    temp = quantity_temp,
    tmax = tmax,
    initial = initial_quantity,
    seed = .magp_design_seed(seed, 1)
  )
  alignment <- .magp_align_initial_design(
    quantity = quantity_search$quantity,
    sequence = sequence_search$sequence,
    weights = joint_weights,
    p = p,
    maxit = alignment_maxit,
    temp = alignment_temp,
    tmax = tmax,
    seed = .magp_design_seed(seed, 2)
  )

  quantity <- alignment$quantity
  sequence <- alignment$sequence
  colnames(quantity) <- paste0("quantity_", seq_len(q))
  colnames(sequence) <- paste0("sequence_", seq_len(q))
  design <- cbind(quantity, sequence)

  output <- list(
    design = design,
    quantity = quantity,
    sequence = sequence,
    criteria = c(
      sequence = sequence_search$criterion,
      quantitative = quantity_search$criterion,
      joint_before_alignment = alignment$initial_criterion,
      joint = alignment$criterion
    ),
    n = n,
    q = q,
    p = p,
    sequence_search = sequence_search,
    sequence_method = sequence_method,
    quantity_search = quantity_search,
    alignment = alignment,
    seed = seed,
    call = call
  )
  class(output) <- "magp_initial_design"
  output
}

#' Print a quantitative-sequence initial design
#'
#' @param x A `magp_initial_design` object.
#' @param ... Additional arguments, currently unused.
#'
#' @return `x`, invisibly.
#'
#' @export
print.magp_initial_design <- function(x, ...) {
  cat("Quantitative-sequence initial design\n")
  cat("  Runs:", x$n, "\n")
  cat("  Components:", x$q, "\n")
  cat("  Sequence method:", x$sequence_search$method, "\n")
  cat("  Sequence criterion:", format(x$criteria[["sequence"]], digits = 7), "\n")
  cat(
    "  Quantitative criterion:",
    format(x$criteria[["quantitative"]], digits = 7), "\n"
  )
  cat("  Joint criterion:", format(x$criteria[["joint"]], digits = 7), "\n")
  invisible(x)
}
