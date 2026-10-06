# Integrated initial-design and Bayesian optimization workflow.

#' Start Bayesian optimization before any experiments have been run
#'
#' Use this function when no initial data are available. It creates an initial
#' quantitative-sequence design, evaluates `FUN` at those runs, and then uses
#' MaGP expected improvement to select later runs.
#'
#' @param FUN Function that evaluates one experiment. It must accept the
#'   generated arguments `quantity_1`, ..., `quantity_q`, `sequence_1`, ...,
#'   `sequence_q`. It must return one finite number, or a list with one finite
#'   number named `Score` or `Value`.
#' @param n_initial Number of experiments in the generated initial design.
#' @param q Number of components. It must be at least three.
#' @param model MaGP mapping to fit. Use `"2d"` for the compact mapping or
#'   `"full"` for the full mapping.
#' @param direction Use `"minimize"` when smaller responses are better and
#'   `"maximize"` when larger responses are better.
#' @param n_iter Maximum number of new experiments selected after the initial
#'   design has been evaluated.
#' @param xi Nonnegative expected-improvement offset. The default, `0`, uses
#'   the current best response as the improvement target.
#' @param stop_ei Nonnegative early-stopping threshold for expected
#'   improvement.
#' @param stop_patience Number of consecutive selected points with expected
#'   improvement less than or equal to `stop_ei` required to stop early.
#' @param seed Optional nonnegative whole-number seed. It makes the initial
#'   design, model starts, acquisition searches, and random values produced by
#'   `FUN` reproducible while preserving the caller's random-number state.
#' @param design_control Optional named list passed to
#'   [magp_initial_design()]. Use `sequence_method` to choose `"random"`,
#'   `"sfta"`, or `"sann"`. Other common choices include `sequence_maxit`,
#'   `sfta_control`, `quantity_maxit`, and `alignment_maxit`. This function
#'   supplies `n`, `q`, and `seed`.
#' @param fit_control Optional named list passed to [magp2d_fit()] or
#'   [magpfull_fit()]. Common choices are `maxeval`, `n_starts`, and `workers`.
#'   This function supplies `X`, `y`, and `seed`.
#' @param acquisition_control Optional named list passed to
#'   [magp_next_point()]. Common choices are `lower`, `upper`, `sequences`,
#'   `n_starts`, `workers`, and `maxit`. This function supplies the fitted
#'   model, direction, `xi`, current best response, and seed.
#' @param objective_args Optional named list of fixed arguments passed to `FUN`
#'   in addition to the generated input columns.
#' @param verbose If `TRUE`, print one line after each initial and sequential
#'   evaluation.
#'
#' @section How to use this function:
#' Define `FUN`, choose `n_initial` and `q`, and state whether the response
#' should be minimized or maximized. The function then:
#'
#' 1. creates an initial design;
#' 2. calls `FUN` once for each initial row;
#' 3. fits the selected MaGP model; and
#' 4. selects and evaluates as many as `n_iter` additional rows.
#'
#' The generated columns are named `quantity_1` through `quantity_q` and
#' `sequence_1` through `sequence_q`. Use [magp_bayes_optimize()] instead when
#' initial experiments and responses already exist.
#'
#' @section Reading the result:
#' The result has the same main components as [magp_bayes_optimize()], including
#' `best_point`, `best_value`, `history`, the final `model`, and all recorded
#' search settings. It also contains:
#'
#' * `initial_design`: the generated quantitative-sequence design;
#' * `initial_response`: the response from each initial run;
#' * `initial_evaluations`: the number of initial runs; and
#' * `design_control`: the initial-design settings that were used; and
#' * `started_from_initial_design`: `TRUE`, marking that the workflow generated
#'   its own starting design.
#'
#' @return An object of class `magp_bayes_opt` containing the initial design,
#'   all completed evaluations, and the final fitted model.
#'
#' @examples
#' \donttest{
#' objective <- function(quantity_1, quantity_2, quantity_3,
#'                       sequence_1, sequence_2, sequence_3) {
#'   quantities <- c(quantity_1, quantity_2, quantity_3)
#'   sequence <- c(sequence_1, sequence_2, sequence_3)
#'   -sum((quantities - c(0.2, 0.6, 0.8))^2) -
#'     0.01 * sum((sequence - c(1, 3, 2))^2)
#' }
#' result <- magp_bayes_optimize_from_scratch(
#'   FUN = objective,
#'   n_initial = 6,
#'   q = 3,
#'   direction = "maximize",
#'   n_iter = 1,
#'   seed = 1,
#'   design_control = list(
#'     sequence_maxit = 100,
#'     quantity_maxit = 100,
#'     alignment_maxit = 100
#'   ),
#'   fit_control = list(maxeval = 100),
#'   acquisition_control = list(n_starts = 2, maxit = 20),
#'   verbose = FALSE
#' )
#' result$initial_design$design
#' result$best_point
#' result$best_value
#' result$history
#' }
#'
#' @export
magp_bayes_optimize_from_scratch <- function(
    FUN, n_initial, q, model = c("2d", "full"),
    direction = c("minimize", "maximize"), n_iter = 10L,
    xi = 0, stop_ei = 0, stop_patience = 3L, seed = NULL,
    design_control = list(), fit_control = list(),
    acquisition_control = list(), objective_args = list(),
    verbose = TRUE) {
  call <- match.call()
  if (!is.function(FUN)) {
    stop("FUN must be a function", call. = FALSE)
  }
  model <- match.arg(model)
  direction <- .magp_match_direction(direction)
  .magp2d_validate_scalar(
    n_initial, "n_initial", lower = 2, integer = TRUE
  )
  .magp2d_validate_scalar(q, "q", lower = 3, integer = TRUE)
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

  design_allowed <- setdiff(
    names(formals(magp_initial_design)),
    c("n", "q", "seed")
  )
  design_control <- .magp_named_control(
    design_control,
    "design_control",
    design_allowed
  )
  fit_function <- if (identical(model, "2d")) magp2d_fit else magpfull_fit
  fit_control <- .magp_named_control(
    fit_control,
    "fit_control",
    setdiff(names(formals(fit_function)), c("X", "y", "seed"))
  )
  acquisition_control <- .magp_named_control(
    acquisition_control,
    "acquisition_control",
    setdiff(
      names(formals(magp_next_point)),
      c("object", "direction", "xi", "best", "seed")
    )
  )
  objective_args <- .magp_named_control(
    objective_args,
    "objective_args",
    names(objective_args)
  )

  input_names <- c(
    paste0("quantity_", seq_len(as.integer(q))),
    paste0("sequence_", seq_len(as.integer(q)))
  )
  duplicated_arguments <- intersect(input_names, names(objective_args))
  if (length(duplicated_arguments)) {
    stop(
      "objective_args duplicates input name(s): ",
      paste(duplicated_arguments, collapse = ", "),
      call. = FALSE
    )
  }

  n_initial <- as.integer(n_initial)
  q <- as.integer(q)
  output <- .magp_with_local_seed(seed, {
    initial_design <- do.call(
      magp_initial_design,
      c(
        list(n = n_initial, q = q, seed = seed),
        design_control
      )
    )

    initial_response <- numeric(n_initial)
    for (i in seq_len(n_initial)) {
      point <- as.data.frame(
        initial_design$design[i, , drop = FALSE],
        check.names = FALSE
      )
      initial_response[[i]] <- .magp_bo_evaluate(
        FUN,
        point,
        objective_args,
        paste0("initial run ", i)
      )
      if (verbose) {
        cat(
          "initial run", i,
          "value", format(initial_response[[i]], digits = 6),
          "\n"
        )
      }
    }

    result <- magp_bayes_optimize(
      FUN = FUN,
      X = initial_design$design,
      y = initial_response,
      model = model,
      direction = direction,
      n_iter = n_iter,
      xi = xi,
      stop_ei = stop_ei,
      stop_patience = stop_patience,
      seed = .magp_seed_offset(seed, 10000L),
      fit_control = fit_control,
      acquisition_control = acquisition_control,
      objective_args = objective_args,
      verbose = verbose
    )
    result$call <- call
    result$initial_design <- initial_design
    result$initial_response <- initial_response
    result$initial_evaluations <- n_initial
    result$design_control <- design_control
    result$started_from_initial_design <- TRUE
    result
  })
  output
}
