# Rcpp-backed fitting and prediction for the two supported mapping models.

.magp_fast_mapping_code <- function(mapping) {
  if (identical(mapping, "2d")) 0L else 1L
}

.magp_fast_npar <- function(q, mapping) {
  if (identical(mapping, "2d")) {
    4L * q - 3L
  } else {
    .magpfull_npar(q)
  }
}

.magp_fast_model_state <- function(X, y, parv, q, tau, mapping) {
  covariance <- cpp_magp_covariance(
    X, parv, q, tau, .magp_fast_mapping_code(mapping)
  )
  cholesky <- chol(covariance)
  inverse_cholesky <- backsolve(cholesky, diag(nrow(X)))
  inverse_covariance <- inverse_cholesky %*% t(inverse_cholesky)
  one <- matrix(1, nrow = nrow(X), ncol = 1L)
  response <- as.matrix(y)
  mean <- as.numeric(
    (t(one) %*% inverse_covariance %*% response) /
      sum(inverse_covariance)
  )
  list(
    mean = mean,
    cholesky = cholesky,
    inverse_covariance = inverse_covariance
  )
}

.magp2d_model_state <- function(X, y, parv, q, tau) {
  .magp_fast_model_state(X, y, parv, q, tau, "2d")
}

.magpfull_model_state <- function(X, y, parv, q, tau) {
  .magp_fast_model_state(X, y, parv, q, tau, "full")
}

.magp_fast_objective_factory <- function(X, y, q, tau, mapping) {
  n <- nrow(X)
  npar <- .magp_fast_npar(q, mapping)
  mapping_code <- .magp_fast_mapping_code(mapping)
  response <- as.matrix(y)
  one <- matrix(1, nrow = n, ncol = 1L)

  function(parv) {
    invalid <- list(objective = 10^10, gradient = rep(0, npar))
    if (length(parv) != npar || any(!is.finite(parv)) ||
        min(parv[seq_len(2L * q)]) <= 0) {
      return(invalid)
    }

    covariance <- cpp_magp_covariance(
      X, parv, q, tau, mapping_code
    )
    cholesky <- try(chol(covariance), silent = TRUE)
    if (inherits(cholesky, "try-error")) {
      return(invalid)
    }
    inverse_cholesky <- backsolve(cholesky, diag(n))
    inverse_covariance <- inverse_cholesky %*% t(inverse_cholesky)

    objective <- 2 * sum(log(diag(cholesky))) +
      t(response) %*% inverse_covariance %*% response -
      (t(one) %*% inverse_covariance %*% response)^2 /
        sum(inverse_covariance)
    mean <- as.numeric(
      t(one) %*% inverse_covariance %*% response /
        sum(inverse_covariance)
    )
    alpha <- as.numeric(
      inverse_covariance %*% (response - mean)
    )
    gradient <- cpp_magp_gradient(
      X, parv, q, mapping_code, inverse_covariance, alpha
    )

    if (!is.finite(objective) || any(!is.finite(gradient))) {
      return(invalid)
    }
    list(
      objective = as.numeric(objective),
      gradient = as.numeric(gradient)
    )
  }
}

.magpfull_objective_factory <- function(X, y, q, tau) {
  .magp_fast_objective_factory(X, y, q, tau, "full")
}

.magp_fast_fit <- function(
    X, y, q, tau, maxeval, xtol_rel,
    lb_sigma, ub_sigma, lb_theta, ub_theta,
    lb_delta, ub_delta, seed, mapping, call) {
  minimum_q <- if (identical(mapping, "2d")) 3 else 2
  prepared <- .magp2d_prepare_fit_inputs(
    X, y, q, minimum_q = minimum_q
  )
  X <- prepared$X
  y <- prepared$y
  q <- prepared$q
  n <- nrow(X)
  npar <- .magp_fast_npar(q, mapping)

  .magp2d_validate_scalar(tau, "tau", lower = 0)
  .magp2d_validate_scalar(maxeval, "maxeval", lower = 1, integer = TRUE)
  .magp2d_validate_scalar(
    xtol_rel, "xtol_rel", lower = 0, strict_lower = TRUE
  )
  .magp2d_validate_scalar(
    lb_sigma, "lb_sigma", lower = 0, strict_lower = TRUE
  )
  .magp2d_validate_scalar(
    ub_sigma, "ub_sigma", lower = 0, strict_lower = TRUE
  )
  .magp2d_validate_scalar(
    lb_theta, "lb_theta", lower = 0, strict_lower = TRUE
  )
  .magp2d_validate_scalar(
    ub_theta, "ub_theta", lower = 0, strict_lower = TRUE
  )
  .magp2d_validate_scalar(lb_delta, "lb_delta")
  .magp2d_validate_scalar(ub_delta, "ub_delta")
  .magp2d_validate_scalar(
    seed, "seed", lower = 0, upper = .Machine$integer.max,
    integer = TRUE, allow_null = TRUE
  )
  if (lb_sigma >= ub_sigma || lb_theta >= ub_theta ||
      lb_delta >= ub_delta) {
    stop(
      "each lower parameter bound must be smaller than its upper bound",
      call. = FALSE
    )
  }
  if (tau == 0 && any(duplicated(as.data.frame(X)))) {
    stop(
      "tau = 0 cannot be used with duplicated training inputs; ",
      "combine duplicates or use a positive nugget variance",
      call. = FALSE
    )
  }

  mapping_count <- npar - 2L * q
  lower <- c(
    rep(lb_sigma, q), rep(lb_theta, q), rep(lb_delta, mapping_count)
  )
  upper <- c(
    rep(ub_sigma, q), rep(ub_theta, q), rep(ub_delta, mapping_count)
  )
  if (!is.null(seed)) set.seed(seed)
  initial <- lower + stats::runif(npar) * (upper - lower)
  evaluate <- .magp_fast_objective_factory(X, y, q, tau, mapping)
  options <- list(
    algorithm = "NLOPT_LD_LBFGS",
    xtol_rel = xtol_rel,
    maxeval = maxeval
  )
  result <- tryCatch(
    nloptr::nloptr(
      x0 = initial,
      eval_f = evaluate,
      lb = lower,
      ub = upper,
      opts = options
    ),
    error = function(e) {
      stop("nloptr failed: ", conditionMessage(e), call. = FALSE)
    }
  )
  valid_result <- is.numeric(result$objective) &&
    length(result$objective) == 1L && is.finite(result$objective) &&
    is.numeric(result$solution) &&
    length(result$solution) == npar &&
    all(is.finite(result$solution))
  if (!valid_result) {
    stop("nloptr returned a non-finite result", call. = FALSE)
  }
  status <- if (is.null(result$status)) NA_integer_ else {
    as.integer(result$status[1L])
  }
  converged <- !is.na(status) && status %in% 1:4
  if (!converged) {
    warning(
      "nloptr did not report convergence; status ", result$status, ": ",
      result$message, ". Inspect the fitted model before using it.",
      call. = FALSE
    )
  }

  solution <- as.vector(result$solution)
  final_state <- .magp_fast_model_state(
    X, y, solution, q, tau, mapping
  )
  parameters <- list(
    sigma2 = solution[seq_len(q)],
    theta = solution[q + seq_len(q)],
    delta = solution[(2L * q + 1L):npar]
  )
  value <- list(
    call = call,
    q = q,
    n = n,
    tau = tau,
    npar = npar,
    X = X,
    y = y,
    input_roles = prepared$input_roles,
    input_names_generated = prepared$input_names_generated,
    quantity_scaling = prepared$quantity_scaling,
    solution = solution,
    mean = final_state$mean,
    parameters = parameters,
    cholesky = final_state$cholesky,
    inverse_covariance = final_state$inverse_covariance,
    objective = as.numeric(result$objective),
    converged = converged,
    initial = initial,
    bounds = list(lower = lower, upper = upper),
    nloptr_status = result,
    engine = "Rcpp"
  )
  if (identical(mapping, "full")) {
    value$mapping <- "full"
    value$mapping_dimension <- q - 1L
    class(value) <- "magpfull"
  } else {
    class(value) <- "magp2d"
  }
  value
}

#' Fit a MaGP model with a two-dimensional sequence map
#'
#' Fits an additive Gaussian process for data that combine component amounts
#' with a component order. Sequence positions are represented by points in a
#' compact two-dimensional latent map.
#'
#' @param X A numeric matrix or data frame. The first `q` input columns contain
#'   quantitative levels and the next `q` columns contain sequence positions.
#'   Each sequence row must be a permutation of `1:q`. A column named `y` may
#'   be included as the response.
#' @param y An optional numeric response vector. It may be omitted when `X`
#'   contains a response column named `y`.
#' @param q The number of components. It is inferred from the number of input
#'   columns when omitted. The two-dimensional model requires at least three
#'   components; the full model requires at least two.
#' @param tau A fixed nonnegative nugget variance added to the covariance
#'   diagonal.
#' @param maxeval Maximum number of objective evaluations used by `nloptr`.
#' @param xtol_rel Relative parameter tolerance used by `nloptr`.
#' @param lb_sigma,ub_sigma Lower and upper bounds for the additive variance
#'   parameters.
#' @param lb_theta,ub_theta Lower and upper bounds for the quantitative
#'   correlation parameters.
#' @param lb_delta,ub_delta Lower and upper bounds for the mapping parameters.
#' @param seed An optional nonnegative integer used to generate the initial
#'   parameter vector.
#'
#' @details The covariance is a sum of component-specific terms. Each term
#'   combines the distance between two quantitative levels with the distance
#'   between their mapped sequence positions. The variance, correlation, and
#'   mapping parameters are estimated together with bounded optimization.
#'
#'   Quantitative columns outside `[0, 1]` are transformed with column-wise
#'   min-max scaling. Their training ranges are stored in the fitted object and
#'   reused for prediction. Columns already in `[0, 1]` are left unchanged.
#'
#' @return An object of class `magp2d`.
#' @examples
#' \donttest{
#' train <- read.table(
#'   system.file("extdata", "example_train.txt", package = "magp"),
#'   header = TRUE
#' )
#' fit <- magp2d_fit(train, seed = 1)
#' fit
#' }
#' @export
magp2d_fit <- function(
    X, y = NULL, q = NULL,
    tau = 0.001, maxeval = 500, xtol_rel = 1e-5,
    lb_sigma = 10, ub_sigma = 1000,
    lb_theta = 0.5, ub_theta = 1000,
    lb_delta = -1, ub_delta = 1,
    seed = NULL) {
  .magp_fast_fit(
    X, y, q, tau, maxeval, xtol_rel,
    lb_sigma, ub_sigma, lb_theta, ub_theta,
    lb_delta, ub_delta, seed, "2d", match.call()
  )
}

#' Fit a MaGP model with a full sequence map
#'
#' Fits the quantitative-sequence model with `q - 1` latent mapping dimensions,
#' giving the sequence positions a less constrained coordinate representation.
#'
#' @inheritParams magp2d_fit
#' @details The data layout, quantitative scaling, and covariance construction
#'   are the same as in [magp2d_fit()]. The difference is the number of latent
#'   coordinates used to represent the sequence positions.
#' @return An object of class `magpfull`.
#' @examples
#' \donttest{
#' train <- read.table(
#'   system.file("extdata", "example_train.txt", package = "magp"),
#'   header = TRUE
#' )
#' fit <- magpfull_fit(train, seed = 1)
#' fit
#' }
#' @export
magpfull_fit <- function(
    X, y = NULL, q = NULL,
    tau = 0.001, maxeval = 500, xtol_rel = 1e-5,
    lb_sigma = 10, ub_sigma = 1000,
    lb_theta = 0.5, ub_theta = 1000,
    lb_delta = -1, ub_delta = 1,
    seed = NULL) {
  .magp_fast_fit(
    X, y, q, tau, maxeval, xtol_rel,
    lb_sigma, ub_sigma, lb_theta, ub_theta,
    lb_delta, ub_delta, seed, "full", match.call()
  )
}

.magp_fast_predict <- function(
    object, newdata, se.fit, type, mapping, expected_class) {
  if (!inherits(object, expected_class)) {
    stop(
      "object must inherit from class '", expected_class, "'",
      call. = FALSE
    )
  }
  if (!is.logical(se.fit) || length(se.fit) != 1L || is.na(se.fit)) {
    stop("se.fit must be TRUE or FALSE", call. = FALSE)
  }
  type <- match.arg(type, c("script", "response", "latent"))

  q <- object$q
  n <- object$n
  tau <- object$tau
  npar <- object$npar
  training <- object$X
  response <- as.matrix(object$y)
  parameters <- object$solution
  expected_npar <- .magp_fast_npar(q, mapping)
  if (!identical(as.integer(npar), as.integer(expected_npar)) ||
      length(parameters) != expected_npar ||
      any(!is.finite(parameters))) {
    stop("object contains an invalid fitted parameter vector", call. = FALSE)
  }
  training_names <- colnames(training)
  prediction_input <- .magp2d_prepare_prediction_inputs(
    newdata,
    q,
    training_names,
    training_names_generated = isTRUE(object$input_names_generated),
    quantity_scaling = object$quantity_scaling
  )

  inverse_covariance <- object$inverse_covariance
  mean <- object$mean
  cholesky <- object$cholesky
  invalid_inverse <- is.null(inverse_covariance) ||
    !identical(dim(inverse_covariance), c(n, n)) ||
    any(!is.finite(inverse_covariance)) ||
    is.null(mean) || length(mean) != 1L || !is.finite(mean)
  invalid_cholesky <- se.fit && (
    is.null(cholesky) || !identical(dim(cholesky), c(n, n)) ||
      any(!is.finite(cholesky))
  )
  if (invalid_inverse || invalid_cholesky) {
    state <- .magp_fast_model_state(
      training, object$y, parameters, q, tau, mapping
    )
    inverse_covariance <- state$inverse_covariance
    mean <- state$mean
    cholesky <- state$cholesky
  }

  cross_covariance <- cpp_magp_cross_covariance(
    training, prediction_input, parameters, q,
    .magp_fast_mapping_code(mapping)
  )
  smooth_self_covariance <- sum(parameters[seq_len(q)])

  if (identical(type, "script")) {
    rounded_training <- round(training, 5L)
    rounded_prediction <- round(prediction_input, 5L)
    for (i in seq_len(nrow(prediction_input))) {
      matches <- rowSums(
        rounded_training != matrix(
          rounded_prediction[i, ],
          nrow = n,
          ncol = 2L * q,
          byrow = TRUE
        )
      ) == 0L
      if (any(matches)) {
        cross_covariance[i, matches] <- smooth_self_covariance + tau
      }
    }
  }

  centered_response <- response - mean
  alpha <- inverse_covariance %*% centered_response
  prediction <- as.numeric(mean + cross_covariance %*% alpha)
  if (!se.fit) return(prediction)

  whitened_one <- forwardsolve(t(cholesky), rep(1, n))
  whitened_gamma <- forwardsolve(t(cholesky), t(cross_covariance))
  quadratic <- colSums(whitened_gamma^2)
  variance_denominator <- sum(whitened_one^2)
  mean_adjustment <- 1 - colSums(whitened_one * whitened_gamma)
  mean_variance <- mean_adjustment^2 / variance_denominator
  target_self_covariance <- smooth_self_covariance +
    if (type %in% c("script", "response")) tau else 0
  raw_variance <- target_self_covariance - quadratic + mean_variance
  tolerance <- sqrt(.Machine$double.eps) * pmax(
    1, abs(target_self_covariance), abs(quadratic), abs(mean_variance)
  )
  invalid <- which(raw_variance < -tolerance)
  if (length(invalid)) {
    i <- invalid[1L]
    stop(
      "materially negative predictive variance at newdata row ", i,
      " under type = '", type,
      "'. This can occur with type = 'script' for duplicated or ",
      "rounded-matching training inputs; use type = 'response' or ",
      "'latent', or inspect the fitted covariance.",
      call. = FALSE
    )
  }
  variance <- pmax(raw_variance, 0)
  list(
    fit = prediction,
    se.fit = sqrt(variance),
    variance = variance,
    type = type
  )
}

#' Predict outcomes from a fitted MaGP model
#'
#' Returns predictions for new quantitative-sequence inputs. Plug-in standard
#' errors and variances can be returned with the predictions.
#'
#' @param object A fitted `magp2d` or `magpfull` object.
#' @param newdata A numeric matrix or data frame with the same quantitative and
#'   sequence inputs used for fitting. A response column named `y` is ignored.
#' @param se.fit Logical; if `TRUE`, return predictive standard errors and
#'   variances in addition to fitted values.
#' @param type Prediction convention. `"script"` uses the fitted training-row
#'   convention when `newdata` exactly matches a training row. `"response"`
#'   includes the nugget variance for a future response, and `"latent"`
#'   returns uncertainty for the noise-free surface.
#' @param ... Additional arguments, currently unused.
#'
#' @return If `se.fit = FALSE`, a numeric vector of predictions. Otherwise, a
#'   list with components `fit`, `se.fit`, `variance`, and `type`.
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
#' predict(fit, test[1:3, ], se.fit = TRUE, type = "response")
#' }
#' @name predict.magp
#' @export
predict.magp2d <- function(
    object, newdata, se.fit = FALSE,
    type = c("script", "response", "latent"), ...) {
  .magp_fast_predict(
    object, newdata, se.fit, type, "2d", "magp2d"
  )
}

#' @rdname predict.magp
#' @export
predict.magpfull <- function(
    object, newdata, se.fit = FALSE,
    type = c("script", "response", "latent"), ...) {
  .magp_fast_predict(
    object, newdata, se.fit, type, "full", "magpfull"
  )
}
