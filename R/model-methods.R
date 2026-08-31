#' Summarize a fitted two-dimensional MaGP model
#'
#' Prints the model size, nugget variance, fitted mean, objective value, and
#' optimizer status.
#'
#' @param x A `magp2d` model returned by [magp2d_fit()].
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.magp2d <- function(x, ...) {
  cat("<magp2d model>\n")
  cat("  components:", x$q, "\n")
  cat("  training runs:", x$n, "\n")
  cat("  nugget:", x$tau, "\n")
  cat("  fitted mean:", round(x$mean, 4), "\n")
  cat("  objective:", round(x$objective, 4), "\n")
  state <- if (isTRUE(x$converged)) "converged" else "not converged"
  cat("  optimizer:", state, "(status", x$nloptr_status$status, ")\n")
  invisible(x)
}

#' Summarize a fitted full-mapping MaGP model
#'
#' Prints the mapping dimension, model size, nugget variance, fitted mean,
#' objective value, and optimizer status.
#'
#' @param x A `magpfull` model returned by [magpfull_fit()].
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.magpfull <- function(x, ...) {
  cat("<magpfull model>\n")
  cat("  mapping dimension:", x$mapping_dimension, "\n")
  cat("  components:", x$q, "\n")
  cat("  training runs:", x$n, "\n")
  cat("  nugget:", x$tau, "\n")
  cat("  fitted mean:", round(x$mean, 4), "\n")
  cat("  objective:", round(x$objective, 4), "\n")
  state <- if (isTRUE(x$converged)) "converged" else "not converged"
  cat("  optimizer:", state, "(status", x$nloptr_status$status, ")\n")
  invisible(x)
}

#' Calculate root-mean-squared prediction error
#'
#' @param pred Numeric vector of predictions.
#' @param actual Numeric vector of observed values with the same length as
#'   `pred`.
#' @description Calculates the root-mean-squared error between predicted and
#'   observed values.
#' @return A single nonnegative numeric value.
#' @export
magp2d_rmse <- function(pred, actual) {
  if (!is.numeric(pred) || is.factor(pred) ||
      !is.numeric(actual) || is.factor(actual)) {
    stop("pred and actual must be numeric vectors", call. = FALSE)
  }
  pred <- as.numeric(pred)
  actual <- as.numeric(actual)
  if (length(pred) != length(actual)) {
    stop("pred and actual must be the same length", call. = FALSE)
  }
  if (!length(pred)) {
    stop("pred and actual must not be empty", call. = FALSE)
  }
  if (any(!is.finite(pred)) || any(!is.finite(actual))) {
    stop("pred and actual must contain only finite values", call. = FALSE)
  }
  sqrt(mean((pred - actual)^2))
}
