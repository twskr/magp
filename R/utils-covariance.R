# Internal covariance and input-processing helpers shared by the two mapping
# models. The covariance calculations use explicit arguments so the same
# implementation can be reused by fitting and prediction.
#
# Here, q is the number of components. For the two-dimensional model, the
# parameter vector contains q variances, q correlation parameters, and
# 2*q - 3 mapping coordinates.

#' Validate a scalar numeric argument
#'
#' @noRd
.magp2d_validate_scalar <- function(x, name, lower = -Inf, upper = Inf,
                                    strict_lower = FALSE, integer = FALSE,
                                    allow_null = FALSE) {
  if (allow_null && is.null(x)) {
    return(invisible(x))
  }
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x)) {
    stop(name, " must be one finite numeric value", call. = FALSE)
  }
  if (integer && x != floor(x)) {
    stop(name, " must be a whole number", call. = FALSE)
  }
  below <- if (strict_lower) x <= lower else x < lower
  if (below || x > upper) {
    interval <- if (strict_lower) "(" else "["
    stop(
      name, " must be in ", interval, lower, ", ", upper, "]",
      call. = FALSE
    )
  }
  invisible(x)
}

#' Validate quantitative-sequence design data
#'
#' @noRd
.magp2d_validate_design <- function(x, q, name, min_rows = 1L,
                                    check_quantity_range = TRUE) {
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop(name, " must be a numeric matrix or data frame", call. = FALSE)
  }
  x <- as.matrix(x)
  if (!is.numeric(x)) {
    stop(name, " must contain only numeric values", call. = FALSE)
  }
  if (ncol(x) != 2L * q) {
    stop(
      name, " must have 2*q input columns: the first q columns must be ",
      "quantitative inputs, and the last q columns must be ",
      "sequence inputs; got ", ncol(x), " columns for q = ", q,
      call. = FALSE
    )
  }
  if (nrow(x) < min_rows) {
    stop(name, " must have at least ", min_rows, " row(s)", call. = FALSE)
  }
  if (any(!is.finite(x))) {
    stop(name, " must contain only finite values (no NA, NaN, or Inf)", call. = FALSE)
  }
  if (!is.null(colnames(x)) && anyDuplicated(colnames(x))) {
    stop(name, " must not have duplicated column names", call. = FALSE)
  }

  names_generated <- is.null(colnames(x)) ||
    any(is.na(colnames(x))) || any(!nzchar(colnames(x)))
  if (names_generated) {
    colnames(x) <- c(
      paste0("quantity_", seq_len(q)),
      paste0("sequence_", seq_len(q))
    )
  }

  quantity_part <- x[, seq_len(q), drop = FALSE]
  invalid_quantity <- which(
    quantity_part < 0 | quantity_part > 1,
    arr.ind = TRUE
  )
  if (isTRUE(check_quantity_range) && nrow(invalid_quantity)) {
    rows <- unique(invalid_quantity[, "row"])
    columns <- unique(invalid_quantity[, "col"])
    shown_rows <- paste(rows[seq_len(min(length(rows), 5L))], collapse = ", ")
    shown_columns <- paste(colnames(quantity_part)[columns], collapse = ", ")
    row_suffix <- if (length(rows) > 5L) ", ..." else ""
    stop(
      name, " must place quantitative inputs in the first half of its ",
      "columns, with every quantitative input scaled to [0, 1]; ",
      "out-of-range value(s) in column(s) ", shown_columns,
      ", row(s) ", shown_rows, row_suffix,
      call. = FALSE
    )
  }

  order_part <- x[, (q + 1L):(2L * q), drop = FALSE]
  valid_order <- apply(order_part, 1L, function(z) {
    all(z == floor(z)) && identical(sort(as.integer(z)), seq_len(q))
  })
  if (any(!valid_order)) {
    bad <- which(!valid_order)
    shown <- paste(bad[seq_len(min(length(bad), 5L))], collapse = ", ")
    suffix <- if (length(bad) > 5L) ", ..." else ""
    stop(
      name, " sequence columns must be a permutation of 1:q in every row; ",
      "invalid row(s): ", shown, suffix,
      call. = FALSE
    )
  }
  attr(x, "magp_names_generated") <- names_generated
  x
}

#' Define the quantitative-input transformation used by a fitted model
#'
#' @noRd
.magp2d_quantity_scaling <- function(x, q) {
  quantity <- x[, seq_len(q), drop = FALSE]
  minimum <- apply(quantity, 2L, min)
  maximum <- apply(quantity, 2L, max)
  span <- maximum - minimum
  tolerance <- sqrt(.Machine$double.eps) * pmax(
    1, abs(minimum), abs(maximum)
  )

  list(
    minimum = as.numeric(minimum),
    maximum = as.numeric(maximum),
    span = as.numeric(span),
    applied = as.logical(minimum < 0 | maximum > 1),
    constant = as.logical(span <= tolerance),
    columns = colnames(quantity),
    method = "auto"
  )
}

#' Apply a fitted quantitative-input transformation
#'
#' @noRd
.magp2d_apply_quantity_scaling <- function(x, q, scaling, name) {
  required <- c(
    "minimum", "maximum", "span", "applied", "constant", "columns"
  )
  if (!is.list(scaling) || !all(required %in% names(scaling)) ||
      any(lengths(scaling[required[1:5]]) != q)) {
    stop("invalid quantitative scaling information", call. = FALSE)
  }

  quantity <- x[, seq_len(q), drop = FALSE]
  for (j in seq_len(q)) {
    values <- quantity[, j]
    tolerance <- sqrt(.Machine$double.eps) * max(
      1, abs(scaling$minimum[j]), abs(scaling$maximum[j])
    )

    if (isTRUE(scaling$constant[j])) {
      bad <- which(abs(values - scaling$minimum[j]) > tolerance)
      if (length(bad)) {
        shown <- paste(bad[seq_len(min(length(bad), 5L))], collapse = ", ")
        suffix <- if (length(bad) > 5L) ", ..." else ""
        stop(
          name, " quantitative column '", colnames(quantity)[j],
          "' was constant in the training data; new values differ at row(s) ",
          shown, suffix,
          call. = FALSE
        )
      }
      if (isTRUE(scaling$applied[j])) {
        quantity[, j] <- 0
      }
      next
    }

    if (isTRUE(scaling$applied[j])) {
      bad <- which(
        values < scaling$minimum[j] - tolerance |
          values > scaling$maximum[j] + tolerance
      )
      if (length(bad)) {
        shown <- paste(bad[seq_len(min(length(bad), 5L))], collapse = ", ")
        suffix <- if (length(bad) > 5L) ", ..." else ""
        stop(
          name, " quantitative column '", colnames(quantity)[j],
          "' must remain within the training range [",
          format(scaling$minimum[j]), ", ", format(scaling$maximum[j]),
          "]; out-of-range row(s): ", shown, suffix,
          call. = FALSE
        )
      }
      quantity[, j] <- (values - scaling$minimum[j]) / scaling$span[j]
    } else {
      bad <- which(values < -tolerance | values > 1 + tolerance)
      if (length(bad)) {
        shown <- paste(bad[seq_len(min(length(bad), 5L))], collapse = ", ")
        suffix <- if (length(bad) > 5L) ", ..." else ""
        stop(
          name, " quantitative column '", colnames(quantity)[j],
          "' must be in [0, 1]; out-of-range row(s): ", shown, suffix,
          call. = FALSE
        )
      }
    }
  }

  quantity <- pmin(pmax(quantity, 0), 1)
  x[, seq_len(q)] <- quantity
  x
}

#' Prepare fitting inputs and infer their roles
#'
#' @noRd
.magp2d_prepare_fit_inputs <- function(X, y = NULL, q = NULL,
                                       minimum_q = 3L) {
  if (!is.matrix(X) && !is.data.frame(X)) {
    stop("X must be a numeric matrix or data frame", call. = FALSE)
  }

  response_columns <- if (is.null(colnames(X))) {
    integer()
  } else {
    which(colnames(X) == "y")
  }
  if (length(response_columns) > 1L) {
    stop("X must contain at most one response column named 'y'", call. = FALSE)
  }
  should_extract_response <- length(response_columns) == 1L && (
    is.null(y) ||
      (is.null(q) && ncol(X) %% 2L == 1L) ||
      (!is.null(q) && ncol(X) == 2L * q + 1L)
  )
  if (should_extract_response) {
    response_from_X <- X[, response_columns, drop = TRUE]
    if (is.null(y)) y <- response_from_X
    X <- X[, -response_columns, drop = FALSE]
  }
  if (is.null(y)) {
    stop(
      "y must be supplied as a numeric vector, or X must contain a ",
      "numeric response column named 'y'",
      call. = FALSE
    )
  }

  if (is.null(q)) {
    if (ncol(X) %% 2L != 0L) {
      stop(
        "X must have an even number of input columns: the first half must ",
        "be quantitative inputs, and the second half must ",
        "be sequence inputs",
        call. = FALSE
      )
    }
    q <- ncol(X) / 2L
  }
  .magp2d_validate_scalar(q, "q", lower = minimum_q, integer = TRUE)
  q <- as.integer(q)
  X <- .magp2d_validate_design(
    X, q, "X", min_rows = 2L, check_quantity_range = FALSE
  )
  names_generated <- isTRUE(attr(X, "magp_names_generated"))
  attr(X, "magp_names_generated") <- NULL
  quantity_scaling <- .magp2d_quantity_scaling(X, q)
  X <- .magp2d_apply_quantity_scaling(X, q, quantity_scaling, "X")

  if (!is.numeric(y) || is.factor(y)) {
    stop("y must be a numeric vector", call. = FALSE)
  }
  y <- as.numeric(y)
  if (length(y) != nrow(X)) {
    stop("length(y) must equal nrow(X)", call. = FALSE)
  }
  if (any(!is.finite(y))) {
    stop(
      "y must contain only finite values (no NA, NaN, or Inf)",
      call. = FALSE
    )
  }

  list(
    X = X,
    y = y,
    q = q,
    input_names_generated = names_generated,
    quantity_scaling = quantity_scaling,
    input_roles = list(
      quantity = colnames(X)[seq_len(q)],
      sequence = colnames(X)[q + seq_len(q)]
    )
  )
}

#' Prepare new quantitative-sequence inputs for prediction
#'
#' @noRd
.magp2d_prepare_prediction_inputs <- function(newdata, q, training_names,
                                              training_names_generated = FALSE,
                                              quantity_scaling = NULL) {
  if (!is.matrix(newdata) && !is.data.frame(newdata)) {
    stop("newdata must be a numeric matrix or data frame", call. = FALSE)
  }

  response_columns <- if (is.null(colnames(newdata))) {
    integer()
  } else {
    which(colnames(newdata) == "y")
  }
  if (length(response_columns) > 1L) {
    stop(
      "newdata must contain at most one response column named 'y'",
      call. = FALSE
    )
  }
  if (length(response_columns) == 1L && ncol(newdata) == 2L * q + 1L) {
    newdata <- newdata[, -response_columns, drop = FALSE]
  }

  prediction_input <- as.matrix(newdata)
  prediction_names <- colnames(prediction_input)
  if (!is.null(training_names) && !is.null(prediction_names)) {
    if (anyDuplicated(prediction_names)) {
      stop("newdata must not have duplicated column names", call. = FALSE)
    }
    if (setequal(training_names, prediction_names)) {
      prediction_input <- prediction_input[, training_names, drop = FALSE]
    } else if (!isTRUE(training_names_generated)) {
      stop(
        "newdata column names must match the training data",
        call. = FALSE
      )
    }
  }
  prediction_input <- .magp2d_validate_design(
    prediction_input, q, "newdata", min_rows = 1L,
    check_quantity_range = is.null(quantity_scaling)
  )
  attr(prediction_input, "magp_names_generated") <- NULL
  if (isTRUE(training_names_generated)) {
    colnames(prediction_input) <- training_names
  }
  if (!is.null(quantity_scaling)) {
    prediction_input <- .magp2d_apply_quantity_scaling(
      prediction_input, q, quantity_scaling, "newdata"
    )
  }
  prediction_input
}
