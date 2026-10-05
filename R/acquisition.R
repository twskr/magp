# Expected improvement and mixed quantitative-sequence acquisition search.

.magp_validate_model <- function(object) {
  if (!inherits(object, c("magp2d", "magpfull"))) {
    stop(
      "object must be a fitted magp2d or magpfull model",
      call. = FALSE
    )
  }
  if (!is.numeric(object$y) || !length(object$y) ||
      any(!is.finite(object$y))) {
    stop("object contains an invalid response vector", call. = FALSE)
  }
  invisible(object)
}

.magp_match_direction <- function(direction) {
  match.arg(direction, c("minimize", "maximize"))
}

.magp_ei_from_prediction <- function(mean, standard_error, best, direction,
                                     xi, variance_floor = 1e-8) {
  if (length(mean) != length(standard_error) ||
      any(!is.finite(mean)) || any(!is.finite(standard_error)) ||
      any(standard_error < 0)) {
    stop(
      "predictive means and standard errors must be finite and compatible",
      call. = FALSE
    )
  }
  improvement <- if (identical(direction, "minimize")) {
    best - mean - xi
  } else {
    mean - best - xi
  }
  usable <- standard_error^2 >= variance_floor
  result <- numeric(length(mean))
  if (any(usable)) {
    z <- improvement[usable] / standard_error[usable]
    result[usable] <-
      improvement[usable] * stats::pnorm(z) +
      standard_error[usable] * stats::dnorm(z)
  }
  pmax(as.numeric(result), 0)
}

#' Score candidate experiments with expected improvement
#'
#' Calculates expected improvement for one or more candidate rows. Higher
#' values indicate candidates that offer a better combination of predicted
#' improvement and uncertainty.
#'
#' @param object A fitted `magp2d` or `magpfull` model.
#' @param newdata A numeric matrix or data frame accepted by the model's
#'   [stats::predict()] method.
#' @param direction Use `"minimize"` when smaller responses are better and
#'   `"maximize"` when larger responses are better.
#' @param best Optional response that a new point should improve upon. When it
#'   is omitted, the function uses the best observed response in `object`.
#' @param xi Nonnegative improvement offset. The default is `0`. Larger values
#'   require a candidate to exceed `best` by more.
#'
#' @details Expected improvement uses the model's latent predictive mean and
#'   standard error. A candidate can receive a high value because its predicted
#'   response is good, its uncertainty is large, or both. Rows already present
#'   in the training data receive a value of zero. Predictions with variance
#'   below `1e-8` also receive zero.
#'
#' @return A nonnegative numeric vector with one expected-improvement value per
#'   row of `newdata`.
#'
#' @references
#' Jones, D. R., Schonlau, M., and Welch, W. J. (1998).
#' Efficient Global Optimization of Expensive Black-Box Functions.
#' Journal of Global Optimization, 13, 455-492.
#' \doi{10.1023/A:1008306431147}.
#'
#' @examples
#' \donttest{
#' train <- read.table(
#'   system.file("extdata", "example_train.txt", package = "magp"),
#'   header = TRUE
#' )
#' test <- read.table(
#'   system.file("extdata", "example_test.txt", package = "magp"),
#'   header = TRUE
#' )
#' fit <- magp2d_fit(train, seed = 1)
#' magp_expected_improvement(
#'   fit,
#'   test[1:3, ],
#'   direction = "minimize"
#' )
#' }
#'
#' @export
magp_expected_improvement <- function(
    object, newdata, direction = c("minimize", "maximize"),
    best = NULL, xi = 0) {
  .magp_validate_model(object)
  direction <- .magp_match_direction(direction)
  .magp2d_validate_scalar(xi, "xi", lower = 0)
  if (is.null(best)) {
    best <- if (identical(direction, "minimize")) {
      min(object$y)
    } else {
      max(object$y)
    }
  }
  .magp2d_validate_scalar(best, "best")

  prediction <- stats::predict(
    object,
    newdata,
    se.fit = TRUE,
    type = "latent"
  )
  improvement <- .magp_ei_from_prediction(
    prediction$fit,
    prediction$se.fit,
    best,
    direction,
    xi
  )
  candidate <- .magp2d_prepare_prediction_inputs(
    newdata,
    object$q,
    colnames(object$X),
    training_names_generated = isTRUE(object$input_names_generated),
    quantity_scaling = object$quantity_scaling
  )
  observed <- vapply(seq_len(nrow(candidate)), function(i) {
    difference <- abs(sweep(
      object$X,
      2L,
      candidate[i, ],
      FUN = "-"
    ))
    any(apply(difference == 0, 1L, all))
  }, logical(1L))
  improvement[observed] <- 0
  improvement
}

.magp_quantity_bounds <- function(object, lower, upper) {
  q <- object$q
  scaling <- object$quantity_scaling
  required <- c("minimum", "maximum", "applied", "constant")
  if (!is.list(scaling) || !all(required %in% names(scaling)) ||
      any(lengths(scaling[required]) != q)) {
    stop("object contains invalid quantitative scaling information",
         call. = FALSE)
  }

  allowed_lower <- ifelse(
    scaling$constant,
    scaling$minimum,
    ifelse(scaling$applied, scaling$minimum, 0)
  )
  allowed_upper <- ifelse(
    scaling$constant,
    scaling$maximum,
    ifelse(scaling$applied, scaling$maximum, 1)
  )

  validate_bound <- function(value, name, default) {
    if (is.null(value)) return(as.numeric(default))
    if (!is.numeric(value) || is.factor(value) ||
        !(length(value) %in% c(1L, q)) || any(!is.finite(value))) {
      stop(name, " must contain one value or q finite numeric values",
           call. = FALSE)
    }
    rep(as.numeric(value), length.out = q)
  }
  lower <- validate_bound(lower, "lower", allowed_lower)
  upper <- validate_bound(upper, "upper", allowed_upper)
  tolerance <- sqrt(.Machine$double.eps) * pmax(
    1, abs(allowed_lower), abs(allowed_upper)
  )

  if (any(lower < allowed_lower - tolerance) ||
      any(upper > allowed_upper + tolerance)) {
    stop(
      "lower and upper must remain within the model's quantitative ",
      "prediction ranges",
      call. = FALSE
    )
  }
  if (any(lower > upper + tolerance)) {
    stop("each lower bound must not exceed its upper bound", call. = FALSE)
  }
  lower <- pmax(lower, allowed_lower)
  upper <- pmin(upper, allowed_upper)
  fixed <- abs(upper - lower) <= tolerance
  upper[fixed] <- lower[fixed]

  quantity_names <- object$input_roles$quantity
  if (is.null(quantity_names) || length(quantity_names) != q) {
    quantity_names <- colnames(object$X)[seq_len(q)]
  }
  names(lower) <- names(upper) <- quantity_names
  list(lower = lower, upper = upper, fixed = fixed)
}

.magp_all_permutations <- function(q) {
  if (q == 1L) return(matrix(1L, nrow = 1L, ncol = 1L))
  previous <- .magp_all_permutations(q - 1L)
  block_size <- nrow(previous)
  result <- matrix(NA_integer_, nrow = q * block_size, ncol = q)
  for (first in seq_len(q)) {
    rows <- (first - 1L) * block_size + seq_len(block_size)
    result[rows, ] <- cbind(
      first,
      previous + (previous >= first)
    )
  }
  result
}

.magp_validate_acquisition_sequences <- function(sequences, q) {
  if (!is.matrix(sequences) && !is.data.frame(sequences)) {
    stop("sequences must be a numeric matrix or data frame", call. = FALSE)
  }
  sequences <- as.matrix(sequences)
  if (!is.numeric(sequences) || ncol(sequences) != q ||
      nrow(sequences) < 1L || any(!is.finite(sequences))) {
    stop(
      "sequences must have q columns and contain finite numeric values",
      call. = FALSE
    )
  }
  valid <- apply(sequences, 1L, function(row) {
    all(row == floor(row)) &&
      identical(sort(as.integer(row)), seq_len(q))
  })
  if (any(!valid)) {
    stop("every row of sequences must be a permutation of 1:q",
         call. = FALSE)
  }
  storage.mode(sequences) <- "integer"
  if (anyDuplicated(as.data.frame(sequences))) {
    stop("sequences must not contain duplicated rows", call. = FALSE)
  }
  colnames(sequences) <- paste0("sequence_", seq_len(q))
  sequences
}

.magp_acquisition_sequences <- function(q, sequences, max_sequences, seed) {
  if (!is.null(sequences)) {
    result <- .magp_validate_acquisition_sequences(sequences, q)
    attr(result, "source") <- "provided"
    return(result)
  }

  if (lgamma(q + 1) <= log(max_sequences) + 1e-12) {
    result <- .magp_all_permutations(q)
    colnames(result) <- paste0("sequence_", seq_len(q))
    attr(result, "source") <- "enumerated"
    return(result)
  }

  result <- .magp_with_local_seed(
    seed,
    .magp_random_sequence_positions(max_sequences, q)
  )
  attr(result, "source") <- "sampled"
  result
}

.magp_seed_offset <- function(seed, offset) {
  if (is.null(seed)) return(NULL)
  as.integer((as.double(seed) + as.double(offset)) %% .Machine$integer.max)
}

.magp_candidate_row <- function(object, quantity, sequence) {
  point <- matrix(
    c(as.numeric(quantity), as.numeric(sequence)),
    nrow = 1L
  )
  colnames(point) <- colnames(object$X)
  point
}

.magp_candidate_is_observed <- function(object, point, tolerance) {
  prepared <- .magp2d_prepare_prediction_inputs(
    point,
    object$q,
    colnames(object$X),
    training_names_generated = isTRUE(object$input_names_generated),
    quantity_scaling = object$quantity_scaling
  )
  difference <- abs(sweep(object$X, 2L, prepared[1L, ], FUN = "-"))
  any(apply(difference <= tolerance, 1L, all))
}

.magp_acquisition_value <- function(object, quantity, sequence, direction,
                                    best, xi, exclude_observed,
                                    duplicate_tolerance) {
  point <- .magp_candidate_row(object, quantity, sequence)
  if (exclude_observed && .magp_candidate_is_observed(
      object, point, duplicate_tolerance)) {
    return(-1)
  }
  magp_expected_improvement(
    object,
    point,
    direction = direction,
    best = best,
    xi = xi
  )[[1L]]
}

.magp_acquisition_tasks <- function(sequences, bounds, n_starts, seed) {
  q <- ncol(sequences)
  midpoint <- (bounds$lower + bounds$upper) / 2
  tasks <- .magp_with_local_seed(seed, {
    result <- vector("list", nrow(sequences) * n_starts)
    task_id <- 0L
    for (sequence_index in seq_len(nrow(sequences))) {
      initial <- matrix(midpoint, nrow = n_starts, ncol = q, byrow = TRUE)
      if (n_starts > 1L) {
        random <- matrix(
          stats::runif((n_starts - 1L) * q),
          nrow = n_starts - 1L,
          ncol = q
        )
        initial[-1L, ] <- matrix(
          bounds$lower,
          nrow = n_starts - 1L,
          ncol = q,
          byrow = TRUE
        ) + random * matrix(
          bounds$upper - bounds$lower,
          nrow = n_starts - 1L,
          ncol = q,
          byrow = TRUE
        )
      }
      for (start in seq_len(n_starts)) {
        task_id <- task_id + 1L
        result[[task_id]] <- list(
          task = task_id,
          sequence_index = sequence_index,
          start = start,
          sequence = sequences[sequence_index, ],
          initial = initial[start, ]
        )
      }
    }
    result
  })
  tasks
}

.magp_optimize_acquisition_task <- function(
    task, object, bounds, direction, best, xi, maxit, factr, pgtol,
    exclude_observed, duplicate_tolerance) {
  variable <- which(!bounds$fixed)
  make_quantity <- function(parameter) {
    quantity <- bounds$lower
    quantity[variable] <- parameter
    quantity
  }
  objective <- function(parameter) {
    value <- .magp_acquisition_value(
      object,
      make_quantity(parameter),
      task$sequence,
      direction,
      best,
      xi,
      exclude_observed,
      duplicate_tolerance
    )
    -value
  }

  result <- tryCatch({
    if (!length(variable)) {
      quantity <- bounds$lower
      value <- .magp_acquisition_value(
        object,
        quantity,
        task$sequence,
        direction,
        best,
        xi,
        exclude_observed,
        duplicate_tolerance
      )
      list(
        quantity = quantity,
        value = value,
        convergence = 0L,
        evaluations = 1L,
        message = "all quantitative inputs were fixed"
      )
    } else {
      optimized <- stats::optim(
        par = task$initial[variable],
        fn = objective,
        method = "L-BFGS-B",
        lower = bounds$lower[variable],
        upper = bounds$upper[variable],
        control = list(maxit = maxit, factr = factr, pgtol = pgtol)
      )
      quantity <- make_quantity(optimized$par)
      value <- .magp_acquisition_value(
        object,
        quantity,
        task$sequence,
        direction,
        best,
        xi,
        exclude_observed,
        duplicate_tolerance
      )
      list(
        quantity = quantity,
        value = value,
        convergence = as.integer(optimized$convergence),
        evaluations = as.integer(optimized$counts[["function"]]),
        message = if (is.null(optimized$message)) "" else optimized$message
      )
    }
  }, error = function(condition) condition)

  if (inherits(result, "error")) {
    return(list(
      task = task$task,
      sequence_index = task$sequence_index,
      start = task$start,
      sequence = task$sequence,
      quantity = rep(NA_real_, object$q),
      ei = NA_real_,
      duplicate = NA,
      convergence = NA_integer_,
      evaluations = NA_integer_,
      message = "",
      error = conditionMessage(result)
    ))
  }

  point <- .magp_candidate_row(object, result$quantity, task$sequence)
  duplicate <- .magp_candidate_is_observed(
    object, point, duplicate_tolerance
  )
  list(
    task = task$task,
    sequence_index = task$sequence_index,
    start = task$start,
    sequence = task$sequence,
    quantity = result$quantity,
    ei = if (result$value < 0) 0 else as.numeric(result$value),
    duplicate = duplicate,
    convergence = result$convergence,
    evaluations = result$evaluations,
    message = result$message,
    error = NA_character_
  )
}

.magp_parallel_acquisition_task <- function(task, arguments) {
  runner <- utils::getFromNamespace(
    ".magp_optimize_acquisition_task",
    "magp"
  )
  do.call(runner, c(list(task = task), arguments))
}
environment(.magp_parallel_acquisition_task) <- baseenv()

.magp_acquisition_diagnostics <- function(results) {
  data.frame(
    task = vapply(results, `[[`, integer(1L), "task"),
    sequence = vapply(results, `[[`, integer(1L), "sequence_index"),
    start = vapply(results, `[[`, integer(1L), "start"),
    expected_improvement = vapply(results, `[[`, numeric(1L), "ei"),
    duplicate = vapply(results, `[[`, logical(1L), "duplicate"),
    convergence = vapply(results, `[[`, integer(1L), "convergence"),
    evaluations = vapply(results, `[[`, integer(1L), "evaluations"),
    message = vapply(results, `[[`, character(1L), "message"),
    error = vapply(results, `[[`, character(1L), "error"),
    stringsAsFactors = FALSE
  )
}

#' Select the next quantitative-sequence experiment
#'
#' Searches the allowed quantitative values and sequence permutations, then
#' returns the unobserved point with the largest expected improvement.
#'
#' @param object A fitted `magp2d` or `magpfull` model.
#' @param direction Use `"minimize"` when smaller responses are better and
#'   `"maximize"` when larger responses are better.
#' @param xi Nonnegative improvement offset used in expected improvement.
#' @param best Optional response that a new point should improve upon. When it
#'   is omitted, the function uses the best observed response in `object`.
#' @param lower,upper Optional lower and upper bounds for the quantitative
#'   inputs. Supply one value for all components or one value per component.
#'   The defaults use the prediction ranges stored in the fitted model.
#' @param sequences Optional matrix of candidate sequence permutations. When
#'   omitted, every permutation is used if their number does not exceed
#'   `max_sequences`; otherwise a reproducible sample is searched.
#' @param max_sequences Largest number of sequence candidates generated when
#'   `sequences` is not supplied.
#' @param n_starts Number of quantitative starting points searched for each
#'   sequence candidate.
#' @param workers Number of local worker processes. Use `1` for sequential
#'   execution. At most two processes are used.
#' @param maxit Maximum optimization iterations for each quantitative start.
#' @param factr,pgtol Advanced convergence settings passed to
#'   [stats::optim()] for its `"L-BFGS-B"` method.
#' @param exclude_observed If `TRUE`, do not return an input that is already in
#'   the training data.
#' @param duplicate_tolerance Nonnegative absolute tolerance used to identify
#'   an observed input after the model's quantitative scaling is applied.
#' @param seed Optional nonnegative whole-number seed for sequence sampling and
#'   quantitative starting points.
#'
#' @section Sequence search:
#' If `sequences` is supplied, only those rows are searched. Otherwise, the
#' function searches every permutation when there are no more than
#' `max_sequences`. For a larger sequence space, it searches a reproducible
#' sample of `max_sequences` permutations when `seed` is supplied.
#'
#' @section Reading the result:
#' The returned object contains:
#'
#' * `call`: the function call;
#' * `point`: the selected input as a one-row data frame;
#' * `expected_improvement`: the score of the selected point;
#' * `predicted_mean` and `predicted_standard_error`: the MaGP prediction;
#' * `direction`, `best_observed`, and `xi`: the expected-improvement settings;
#' * `bounds`: the quantitative lower and upper bounds;
#' * `sequences` and `sequence_source`: the permutations searched and how they
#'   were obtained;
#' * `selected_sequence` and `selected_start`: the winning search indices;
#' * `n_starts`: the number of quantitative starts per sequence;
#' * `workers_requested`, `workers_used`, and `execution`: the parallel-search
#'   settings actually used; and
#' * `diagnostics`: the result of every sequence and starting-point search.
#'
#' @return An object of class `magp_next_point` containing the selected point,
#'   its prediction, and search diagnostics.
#'
#' @examples
#' \donttest{
#' train <- read.table(
#'   system.file("extdata", "example_train.txt", package = "magp"),
#'   header = TRUE
#' )
#' fit <- magp2d_fit(train, seed = 1)
#' next_run <- magp_next_point(
#'   fit,
#'   direction = "minimize",
#'   n_starts = 2,
#'   maxit = 20,
#'   seed = 2
#' )
#' next_run$point
#' next_run$expected_improvement
#' }
#'
#' @export
magp_next_point <- function(
    object, direction = c("minimize", "maximize"), xi = 0,
    best = NULL, lower = NULL, upper = NULL, sequences = NULL,
    max_sequences = 120L, n_starts = 5L, workers = 1L,
    maxit = 100L, factr = 1e7, pgtol = 0,
    exclude_observed = TRUE,
    duplicate_tolerance = sqrt(.Machine$double.eps), seed = NULL) {
  call <- match.call()
  .magp_validate_model(object)
  direction <- .magp_match_direction(direction)
  .magp2d_validate_scalar(xi, "xi", lower = 0)
  .magp2d_validate_scalar(
    max_sequences, "max_sequences", lower = 1, integer = TRUE
  )
  .magp2d_validate_scalar(n_starts, "n_starts", lower = 1, integer = TRUE)
  .magp2d_validate_scalar(workers, "workers", lower = 1, integer = TRUE)
  .magp2d_validate_scalar(maxit, "maxit", lower = 1, integer = TRUE)
  .magp2d_validate_scalar(factr, "factr", lower = 0, strict_lower = TRUE)
  .magp2d_validate_scalar(pgtol, "pgtol", lower = 0)
  .magp2d_validate_scalar(
    duplicate_tolerance,
    "duplicate_tolerance",
    lower = 0
  )
  .magp2d_validate_scalar(
    seed, "seed", lower = 0, upper = .Machine$integer.max,
    integer = TRUE, allow_null = TRUE
  )
  if (!is.logical(exclude_observed) || length(exclude_observed) != 1L ||
      is.na(exclude_observed)) {
    stop("exclude_observed must be TRUE or FALSE", call. = FALSE)
  }
  if (is.null(best)) {
    best <- if (identical(direction, "minimize")) {
      min(object$y)
    } else {
      max(object$y)
    }
  }
  .magp2d_validate_scalar(best, "best")

  max_sequences <- as.integer(max_sequences)
  n_starts <- as.integer(n_starts)
  requested_workers <- as.integer(workers)
  maxit <- as.integer(maxit)
  bounds <- .magp_quantity_bounds(object, lower, upper)
  candidate_sequences <- .magp_acquisition_sequences(
    object$q,
    sequences,
    max_sequences,
    .magp_seed_offset(seed, 1L)
  )
  sequence_source <- attr(candidate_sequences, "source")
  attr(candidate_sequences, "source") <- NULL
  tasks <- .magp_acquisition_tasks(
    candidate_sequences,
    bounds,
    n_starts,
    .magp_seed_offset(seed, 2L)
  )
  arguments <- list(
    object = object,
    bounds = bounds,
    direction = direction,
    best = best,
    xi = xi,
    maxit = maxit,
    factr = factr,
    pgtol = pgtol,
    exclude_observed = exclude_observed,
    duplicate_tolerance = duplicate_tolerance
  )
  workers_used <- min(requested_workers, length(tasks), 2L)

  if (workers_used == 1L) {
    results <- lapply(tasks, function(task) {
      do.call(
        .magp_optimize_acquisition_task,
        c(list(task = task), arguments)
      )
    })
    execution <- "sequential"
  } else {
    cluster <- tryCatch(
      parallel::makePSOCKcluster(workers_used),
      error = function(condition) {
        stop(
          "could not start parallel workers: ",
          conditionMessage(condition),
          call. = FALSE
        )
      }
    )
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    library_paths <- .libPaths()
    parallel::clusterCall(cluster, function(paths) {
      .libPaths(paths)
      loadNamespace("magp")
      NULL
    }, library_paths)
    results <- parallel::parLapply(
      cluster,
      tasks,
      .magp_parallel_acquisition_task,
      arguments = arguments
    )
    execution <- "PSOCK"
  }

  diagnostics <- .magp_acquisition_diagnostics(results)
  eligible <- which(
    is.finite(diagnostics$expected_improvement) &
      (!exclude_observed | !diagnostics$duplicate)
  )
  if (!length(eligible)) {
    reasons <- unique(diagnostics$error[
      !is.na(diagnostics$error) & nzchar(diagnostics$error)
    ])
    stop(
      "the acquisition search did not produce an eligible point",
      if (length(reasons)) paste0(": ", paste(reasons, collapse = "; ")),
      call. = FALSE
    )
  }
  best_task <- eligible[
    which.max(diagnostics$expected_improvement[eligible])
  ]
  selected <- results[[best_task]]
  point <- .magp_candidate_row(
    object,
    selected$quantity,
    selected$sequence
  )
  point <- as.data.frame(point, check.names = FALSE)
  prediction <- stats::predict(
    object,
    point,
    se.fit = TRUE,
    type = "latent"
  )
  expected_improvement <- magp_expected_improvement(
    object,
    point,
    direction = direction,
    best = best,
    xi = xi
  )[[1L]]

  value <- list(
    call = call,
    point = point,
    expected_improvement = expected_improvement,
    predicted_mean = prediction$fit[[1L]],
    predicted_standard_error = prediction$se.fit[[1L]],
    direction = direction,
    best_observed = best,
    xi = xi,
    bounds = bounds[c("lower", "upper")],
    sequences = candidate_sequences,
    sequence_source = sequence_source,
    selected_sequence = selected$sequence_index,
    selected_start = selected$start,
    n_starts = n_starts,
    workers_requested = requested_workers,
    workers_used = workers_used,
    execution = execution,
    diagnostics = diagnostics
  )
  class(value) <- "magp_next_point"
  value
}

#' Print a MaGP acquisition-search result
#'
#' @param x A `magp_next_point` object returned by [magp_next_point()].
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.magp_next_point <- function(x, ...) {
  cat("<magp next point>\n")
  cat("  direction:", x$direction, "\n")
  cat("  expected improvement:",
      format(x$expected_improvement, digits = 6), "\n")
  cat("  predicted mean:", format(x$predicted_mean, digits = 6), "\n")
  cat("  predicted standard error:",
      format(x$predicted_standard_error, digits = 6), "\n")
  cat("  sequences searched:", nrow(x$sequences),
      paste0("(", x$sequence_source, ")"), "\n")
  cat("  execution:", x$execution, "with", x$workers_used,
      "worker(s)\n")
  invisible(x)
}
