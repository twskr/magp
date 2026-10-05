# Sequential Bayesian optimization with a MaGP surrogate.

.magp_named_control <- function(control, name, allowed) {
  if (!is.list(control)) {
    stop(name, " must be a list", call. = FALSE)
  }
  if (!length(control)) return(control)
  control_names <- names(control)
  if (is.null(control_names) || any(!nzchar(control_names)) ||
      anyDuplicated(control_names)) {
    stop(name, " must have unique, nonempty names", call. = FALSE)
  }
  unknown <- setdiff(control_names, allowed)
  if (length(unknown)) {
    stop(
      name, " contains unsupported argument(s): ",
      paste(unknown, collapse = ", "),
      call. = FALSE
    )
  }
  control
}

.magp_bo_prepare_initial_data <- function(X, y) {
  if (!is.matrix(X) && !is.data.frame(X)) {
    stop("X must be a numeric matrix or data frame", call. = FALSE)
  }
  response_columns <- if (is.null(colnames(X))) {
    integer()
  } else {
    which(colnames(X) == "y")
  }
  if (length(response_columns) > 1L) {
    stop("X must contain at most one response column named 'y'",
         call. = FALSE)
  }
  if (length(response_columns) == 1L) {
    embedded_response <- X[, response_columns, drop = TRUE]
    X <- X[, -response_columns, drop = FALSE]
    if (is.null(y)) {
      y <- embedded_response
    } else if (!isTRUE(all.equal(
      as.numeric(y),
      as.numeric(embedded_response),
      check.attributes = FALSE
    ))) {
      stop("y does not match the response column in X", call. = FALSE)
    }
  }
  if (is.null(y)) {
    stop(
      "y must be supplied, or X must contain a response column named 'y'",
      call. = FALSE
    )
  }
  if (ncol(X) %% 2L != 0L) {
    stop("X must contain an even number of input columns", call. = FALSE)
  }
  q <- ncol(X) / 2L
  if (q < 2L) {
    stop("X must describe at least two components", call. = FALSE)
  }
  X <- .magp2d_validate_design(
    X,
    as.integer(q),
    "X",
    min_rows = 2L,
    check_quantity_range = FALSE
  )
  attr(X, "magp_names_generated") <- NULL
  if (!is.numeric(y) || is.factor(y)) {
    stop("y must be a numeric vector", call. = FALSE)
  }
  y <- as.numeric(y)
  if (length(y) != nrow(X) || any(!is.finite(y))) {
    stop(
      "y must contain one finite numeric value for every row of X",
      call. = FALSE
    )
  }
  list(X = X, y = y, q = as.integer(q))
}

.magp_bo_fit <- function(model, X, y, control, seed) {
  fit_function <- if (identical(model, "2d")) {
    magp2d_fit
  } else {
    magpfull_fit
  }
  arguments <- c(list(X = X, y = y, seed = seed), control)
  do.call(fit_function, arguments)
}

.magp_bo_evaluate <- function(FUN, point, objective_args, context) {
  point_arguments <- as.list(as.numeric(point[1L, ]))
  names(point_arguments) <- colnames(point)
  duplicated_arguments <- intersect(names(point_arguments), names(objective_args))
  if (length(duplicated_arguments)) {
    stop(
      "objective_args duplicates input name(s): ",
      paste(duplicated_arguments, collapse = ", "),
      call. = FALSE
    )
  }
  result <- tryCatch(
    do.call(FUN, c(point_arguments, objective_args)),
    error = function(condition) condition
  )
  if (inherits(result, "error")) {
    stop(
      "FUN failed at ", context, ": ",
      conditionMessage(result),
      call. = FALSE
    )
  }
  value <- if (is.numeric(result) && !is.factor(result) &&
               length(result) == 1L) {
    result
  } else if (is.list(result) && "Score" %in% names(result)) {
    result$Score
  } else if (is.list(result) && "Value" %in% names(result)) {
    result$Value
  } else {
    stop(
      "FUN must return one numeric value or a list containing Score or Value",
      call. = FALSE
    )
  }
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
    stop("FUN returned a non-finite or non-scalar value", call. = FALSE)
  }
  as.numeric(value)
}

.magp_bo_initial_history <- function(X, y) {
  data.frame(
    Iteration = rep(0L, nrow(X)),
    Initial = rep(TRUE, nrow(X)),
    as.data.frame(X, check.names = FALSE),
    Value = y,
    ExpectedImprovement = rep(NA_real_, nrow(X)),
    PredictedMean = rep(NA_real_, nrow(X)),
    PredictedStandardError = rep(NA_real_, nrow(X)),
    check.names = FALSE
  )
}

.magp_bo_history_row <- function(iteration, point, value, acquisition) {
  data.frame(
    Iteration = as.integer(iteration),
    Initial = FALSE,
    point,
    Value = value,
    ExpectedImprovement = acquisition$expected_improvement,
    PredictedMean = acquisition$predicted_mean,
    PredictedStandardError = acquisition$predicted_standard_error,
    check.names = FALSE
  )
}

#' Continue Bayesian optimization from completed experiments
#'
#' Use this function when initial experiments and their responses are already
#' available. It fits a MaGP model, selects a new input with expected
#' improvement, evaluates `FUN`, and adds the new result to the data. This
#' process repeats for at most `n_iter` new evaluations.
#'
#' @param FUN Function that evaluates one experiment. Its argument names must
#'   match the columns of `X`. It must return one finite number, or a list with
#'   one finite number named `Score` or `Value`.
#' @param X Initial inputs as a numeric matrix or data frame. The first half of
#'   the columns contains quantitative values. The second half contains the
#'   sequence positions, with each row forming a permutation of `1:q`. `X` may
#'   also contain a response column named `y`.
#' @param y Numeric responses for the rows of `X`. Leave this as `NULL` when
#'   `X` contains a column named `y`.
#' @param model MaGP mapping to fit. Use `"2d"` for the compact mapping or
#'   `"full"` for the full mapping.
#' @param direction Use `"minimize"` when smaller responses are better and
#'   `"maximize"` when larger responses are better.
#' @param n_iter Maximum number of new experiments to evaluate.
#' @param xi Nonnegative expected-improvement offset. The default, `0`, uses
#'   the current best response as the improvement target. Larger values require
#'   a candidate to exceed that target by more.
#' @param stop_ei Nonnegative early-stopping threshold for expected
#'   improvement.
#' @param stop_patience Number of consecutive selected points with expected
#'   improvement less than or equal to `stop_ei` required to stop early.
#' @param seed Optional nonnegative whole-number seed for reproducible model
#'   starts and acquisition searches. The caller's random-number state is
#'   preserved.
#' @param fit_control Optional named list passed to [magp2d_fit()] or
#'   [magpfull_fit()]. Common choices include `tau`, `maxeval`, `n_starts`, and
#'   `workers`. This function supplies `X`, `y`, and `seed`.
#' @param acquisition_control Optional named list passed to
#'   [magp_next_point()]. Common choices include `lower`, `upper`, `sequences`,
#'   `n_starts`, `workers`, and `maxit`. This function supplies the fitted
#'   model, direction, `xi`, current best response, and seed.
#' @param objective_args Optional named list of fixed arguments passed to `FUN`
#'   in addition to the input columns.
#' @param verbose If `TRUE`, print the observed response and expected
#'   improvement after each new evaluation.
#'
#' @section How the loop works:
#' The initial rows of `X` are treated as completed experiments and are not
#' evaluated again. Each iteration performs four steps:
#'
#' 1. fit the selected MaGP model to all results collected so far;
#' 2. use [magp_next_point()] to select an unobserved input;
#' 3. call `FUN` at that input; and
#' 4. add the response and refit the model.
#'
#' `FUN` receives one named argument for each input column. For columns `A`,
#' `B`, `a`, and `b`, for example, the call is equivalent to
#' `FUN(A = ..., B = ..., a = ..., b = ...)`.
#'
#' @section Reading the result:
#' The returned object contains:
#'
#' * `call`: the function call;
#' * `best_point`: the input row with the best observed response;
#' * `best_value`: the best observed response;
#' * `best_index`: the row number of the best result in `X` and `y`;
#' * `history`: the initial and newly evaluated rows in evaluation order;
#' * `model`: the final fitted MaGP model;
#' * `X` and `y`: all inputs and responses used by the final model;
#' * `acquisitions`: details from each call to [magp_next_point()];
#' * `mapping` and `direction`: the model and optimization direction;
#' * `iterations_requested` and `iterations_completed`: the requested and
#'   completed numbers of new evaluations;
#' * `stop_reason`: why the loop ended;
#' * `xi`, `stop_ei`, and `stop_patience`: the acquisition and stopping
#'   settings; and
#' * `fit_control` and `acquisition_control`: the control lists used in the
#'   search.
#'
#' @return An object of class `magp_bayes_opt`. The final fitted model includes
#'   every completed evaluation.
#'
#' @references
#' Jones, D. R., Schonlau, M., and Welch, W. J. (1998).
#' Efficient Global Optimization of Expensive Black-Box Functions.
#' Journal of Global Optimization, 13, 455-492.
#' \doi{10.1023/A:1008306431147}.
#'
#' @examples
#' \donttest{
#' design <- magp_initial_design(n = 6, q = 3, seed = 1)
#' objective <- function(quantity_1, quantity_2, quantity_3,
#'                       sequence_1, sequence_2, sequence_3) {
#'   quantities <- c(quantity_1, quantity_2, quantity_3)
#'   positions <- c(sequence_1, sequence_2, sequence_3)
#'   -sum((quantities - c(0.2, 0.6, 0.8))^2) -
#'     0.02 * sum((positions - c(1, 3, 2))^2)
#' }
#' initial_y <- apply(design$design, 1L, function(row) {
#'   do.call(objective, as.list(row))
#' })
#' result <- magp_bayes_optimize(
#'   FUN = objective,
#'   X = design$design,
#'   y = initial_y,
#'   direction = "maximize",
#'   n_iter = 1,
#'   seed = 2,
#'   fit_control = list(maxeval = 100),
#'   acquisition_control = list(n_starts = 2, maxit = 20),
#'   verbose = FALSE
#' )
#' result$best_point
#' result$best_value
#' result$history
#' }
#'
#' @export
magp_bayes_optimize <- function(
    FUN, X, y = NULL, model = c("2d", "full"),
    direction = c("minimize", "maximize"), n_iter = 10L,
    xi = 0, stop_ei = 0, stop_patience = 3L, seed = NULL,
    fit_control = list(), acquisition_control = list(),
    objective_args = list(), verbose = TRUE) {
  call <- match.call()
  if (!is.function(FUN)) {
    stop("FUN must be a function", call. = FALSE)
  }
  model <- match.arg(model)
  direction <- .magp_match_direction(direction)
  .magp2d_validate_scalar(n_iter, "n_iter", lower = 0, integer = TRUE)
  .magp2d_validate_scalar(xi, "xi", lower = 0)
  .magp2d_validate_scalar(stop_ei, "stop_ei", lower = 0)
  .magp2d_validate_scalar(
    stop_patience, "stop_patience", lower = 1, integer = TRUE
  )
  .magp2d_validate_scalar(
    seed, "seed", lower = 0, upper = .Machine$integer.max,
    integer = TRUE, allow_null = TRUE
  )
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be TRUE or FALSE", call. = FALSE)
  }

  fit_function <- if (identical(model, "2d")) magp2d_fit else magpfull_fit
  fit_allowed <- setdiff(names(formals(fit_function)), c("X", "y", "seed"))
  fit_control <- .magp_named_control(
    fit_control,
    "fit_control",
    fit_allowed
  )
  acquisition_allowed <- setdiff(
    names(formals(magp_next_point)),
    c("object", "direction", "xi", "best", "seed")
  )
  acquisition_control <- .magp_named_control(
    acquisition_control,
    "acquisition_control",
    acquisition_allowed
  )
  objective_args <- .magp_named_control(
    objective_args,
    "objective_args",
    names(objective_args)
  )
  data <- .magp_bo_prepare_initial_data(X, y)
  if (identical(model, "2d") && data$q < 3L) {
    stop("model = '2d' requires at least three components", call. = FALSE)
  }
  if (length(intersect(colnames(data$X), names(objective_args)))) {
    stop("objective_args must not duplicate an input column name",
         call. = FALSE)
  }

  n_iter <- as.integer(n_iter)
  stop_patience <- as.integer(stop_patience)
  output <- .magp_with_local_seed(seed, {
    current_X <- data$X
    current_y <- data$y
    fit <- .magp_bo_fit(
      model,
      current_X,
      current_y,
      fit_control,
      .magp_seed_offset(seed, 100L)
    )
    history <- .magp_bo_initial_history(current_X, current_y)
    acquisitions <- vector("list", n_iter)
    completed <- 0L
    low_ei_count <- 0L
    stop_reason <- "maximum iterations reached"

    if (n_iter > 0L) {
      for (iteration in seq_len(n_iter)) {
        reference <- if (identical(direction, "minimize")) {
          min(current_y)
        } else {
          max(current_y)
        }
        acquisition <- do.call(
          magp_next_point,
          c(
            list(
              object = fit,
              direction = direction,
              xi = xi,
              best = reference,
              seed = .magp_seed_offset(seed, 1000L + iteration)
            ),
            acquisition_control
          )
        )
        value <- .magp_bo_evaluate(
          FUN,
          acquisition$point,
          objective_args,
          paste0("iteration ", iteration)
        )
        current_X <- rbind(current_X, as.matrix(acquisition$point))
        current_y <- c(current_y, value)
        history <- rbind(
          history,
          .magp_bo_history_row(iteration, acquisition$point, value, acquisition)
        )
        acquisitions[[iteration]] <- acquisition
        completed <- iteration

        fit <- .magp_bo_fit(
          model,
          current_X,
          current_y,
          fit_control,
          .magp_seed_offset(seed, 2000L + iteration)
        )
        if (verbose) {
          cat(
            "iteration", iteration,
            "value", format(value, digits = 6),
            "EI", format(acquisition$expected_improvement, digits = 6),
            "\n"
          )
        }

        if (acquisition$expected_improvement <= stop_ei) {
          low_ei_count <- low_ei_count + 1L
        } else {
          low_ei_count <- 0L
        }
        if (low_ei_count >= stop_patience) {
          stop_reason <- "expected improvement threshold reached"
          break
        }
      }
    }
    acquisitions <- acquisitions[seq_len(completed)]
    best_index <- if (identical(direction, "minimize")) {
      which.min(current_y)
    } else {
      which.max(current_y)
    }
    best_point <- as.numeric(current_X[best_index, ])
    names(best_point) <- colnames(current_X)

    value <- list(
      call = call,
      best_point = best_point,
      best_value = current_y[[best_index]],
      best_index = best_index,
      history = history,
      model = fit,
      X = current_X,
      y = current_y,
      acquisitions = acquisitions,
      mapping = model,
      direction = direction,
      iterations_requested = n_iter,
      iterations_completed = completed,
      stop_reason = stop_reason,
      xi = xi,
      stop_ei = stop_ei,
      stop_patience = stop_patience,
      fit_control = fit_control,
      acquisition_control = acquisition_control
    )
    class(value) <- "magp_bayes_opt"
    value
  })
  output
}

#' Print a MaGP Bayesian optimization result
#'
#' @param x A `magp_bayes_opt` object returned by
#'   [magp_bayes_optimize()].
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.magp_bayes_opt <- function(x, ...) {
  cat("<magp Bayesian optimization>\n")
  cat("  mapping:", x$mapping, "\n")
  cat("  direction:", x$direction, "\n")
  if (isTRUE(x$started_from_initial_design)) {
    cat("  initial design runs:", x$initial_evaluations, "\n")
  }
  cat("  new evaluations:", x$iterations_completed, "\n")
  cat("  best value:", format(x$best_value, digits = 6), "\n")
  cat("  stopped:", x$stop_reason, "\n")
  invisible(x)
}
