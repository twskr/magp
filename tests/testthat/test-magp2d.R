example_data <- function() {
  train <- read.table(
    system.file("extdata", "example_train.txt", package = "magp"),
    header = TRUE
  )
  test <- read.table(
    system.file("extdata", "example_test.txt", package = "magp"),
    header = TRUE
  )
  cols <- c("A", "B", "C", "D", "a", "b", "c", "d")
  list(
    train = train,
    test = test,
    X = as.matrix(train[, cols]),
    y = train$y,
    X_test = as.matrix(test[, cols]),
    y_test = test$y,
    cols = cols
  )
}

reference_solution <- c(
  242.1324, 335.6161, 536.0240, 873.3480,
  193.0639, 897.5149, 942.4635, 652.4657,
  -3.015898e-03, -7.592647e-04, 1.372593e-04,
  -5.599387e-03, -9.738390e-04
)

reference_model <- function(dat, tau = 0.001) {
  structure(
    list(
      q = 4L,
      n = nrow(dat$X),
      tau = tau,
      npar = 13L,
      X = dat$X,
      y = dat$y,
      solution = reference_solution
    ),
    class = "magp2d"
  )
}

uncertainty_fixture <- function(tau = 0.05) {
  X <- rbind(
    c(0, 0, 0, 1, 2, 3),
    c(1, 0, 0, 1, 2, 3),
    c(0, 1, 0, 2, 1, 3),
    c(0, 0, 1, 3, 2, 1)
  )
  colnames(X) <- c("x1", "x2", "x3", "o1", "o2", "o3")
  structure(
    list(
      q = 3L,
      n = 4L,
      tau = tau,
      npar = 9L,
      X = X,
      y = c(2, -1, 4, 0.5),
      solution = c(1.3, 0.7, 0.4, 1.1, 0.9, 0.8, 0.2, -0.3, 0.4)
    ),
    class = "magp2d"
  )
}

test_that("two-dimensional analytical gradients match finite differences", {
  dat <- example_data()
  X <- dat$X[1:8, , drop = FALSE]
  y <- dat$y[1:8]
  evaluate <- magp:::.magp_fast_objective_factory(
    X, y, 4L, 0.001, "2d"
  )
  parameter <- c(
    120, 90, 110, 80,
    4, 7, 3, 5,
    0.2, -0.3, 0.4, -0.1, 0.25
  )
  analytical <- evaluate(parameter)$gradient
  step <- 1e-6
  numerical <- vapply(seq_along(parameter), function(i) {
    plus <- parameter
    minus <- parameter
    plus[i] <- plus[i] + step
    minus[i] <- minus[i] - step
    (
      evaluate(plus)$objective - evaluate(minus)$objective
    ) / (2 * step)
  }, numeric(1))
  scaled_error <- abs(analytical - numerical) /
    pmax(1, abs(analytical), abs(numerical))

  expect_length(analytical, 13L)
  expect_lt(max(scaled_error), 1e-6)
})

test_that("fit and prediction reproduce the stored numerical reference", {
  dat <- example_data()
  fit <- magp2d_fit(dat$X, dat$y, q = 4, tau = 0.001, seed = 1)

  expect_s3_class(fit, "magp2d")
  expect_true(fit$converged)
  expect_equal(fit$q, 4L)
  expect_equal(fit$n, nrow(dat$X))
  expect_length(fit$solution, 13L)
  expect_equal(fit$solution, reference_solution, tolerance = 1e-4)
  expect_true(is.finite(fit$mean))
  expect_named(fit$parameters, c("sigma2", "theta", "delta"))
  expect_length(fit$parameters$sigma2, 4L)
  expect_length(fit$parameters$theta, 4L)
  expect_length(fit$parameters$delta, 5L)
  expect_length(fit$initial, 13L)

  pred <- predict(fit, dat$X_test)
  expect_identical(
    pred,
    predict(fit, dat$X_test, se.fit = FALSE, type = "script")
  )
  expect_length(pred, nrow(dat$X_test))
  expect_true(all(is.finite(pred)))
  expect_equal(magp2d_rmse(pred, dat$y_test), 11.44883, tolerance = 1e-4)
})

test_that("fit validates model inputs before fitting", {
  dat <- example_data()
  X <- dat$X[1:6, , drop = FALSE]
  y <- dat$y[1:6]

  expect_error(magp2d_fit(X, y, q = 2), "q")
  expect_error(magp2d_fit(as.character(X), y, q = 4), "numeric matrix")
  expect_error(magp2d_fit(X, factor(y), q = 4), "numeric vector")

  X_na <- X
  X_na[1, 1] <- NA_real_
  expect_error(magp2d_fit(X_na, y, q = 4), "finite")

  y_na <- y
  y_na[1] <- NA_real_
  expect_error(magp2d_fit(X, y_na, q = 4), "finite")

  X_bad_order <- X
  X_bad_order[1, 5:8] <- c(1, 1, 3, 4)
  expect_error(magp2d_fit(X_bad_order, y, q = 4), "permutation")

  expect_error(magp2d_fit(X, y, q = 4, tau = -0.1), "tau")
  expect_error(magp2d_fit(X, y, q = 4, maxeval = 0), "maxeval")
  expect_error(magp2d_fit(X, y, q = 4, xtol_rel = 0), "xtol_rel")
  expect_error(
    magp2d_fit(X, y, q = 4, lb_sigma = 100, ub_sigma = 10),
    "lower parameter bound"
  )
  expect_error(magp2d_fit(X, y, q = 4, seed = -1), "seed")
})

test_that("fit infers response, q, and input roles without manual columns", {
  dat <- example_data()
  expect_warning(
    fit <- magp2d_fit(dat$train, seed = 1, maxeval = 1),
    "did not report convergence"
  )

  expect_equal(fit$q, 4L)
  expect_equal(fit$y, dat$train$y)
  expect_equal(colnames(fit$X), dat$cols)
  expect_equal(fit$input_roles$quantity, c("A", "B", "C", "D"))
  expect_equal(fit$input_roles$sequence, c("a", "b", "c", "d"))
  expect_named(fit$input_roles, c("quantity", "sequence"))
  expect_false(any(fit$quantity_scaling$applied))
  expect_false(fit$input_names_generated)
  expect_equal(
    predict(fit, dat$test[1:3, ]),
    predict(fit, dat$X_test[1:3, , drop = FALSE])
  )

  unnamed <- dat$X[1:6, , drop = FALSE]
  colnames(unnamed) <- NULL
  expect_warning(
    unnamed_fit <- magp2d_fit(
      unnamed, dat$y[1:6], seed = 1, maxeval = 1
    ),
    "did not report convergence"
  )
  expect_true(unnamed_fit$input_names_generated)
  expect_equal(
    colnames(unnamed_fit$X),
    c(
      paste0("quantity_", 1:4),
      paste0("sequence_", 1:4)
    )
  )
})

test_that("quantitative inputs are scaled once and reused for prediction", {
  dat <- example_data()
  X <- dat$X[1:6, , drop = FALSE]
  y <- dat$y[1:6]

  normalized <- X
  normalized[, 1:4] <- cbind(
    seq(0, 1, length.out = 6),
    c(0, 0.2, 0.4, 0.6, 0.8, 1),
    c(1, 0.8, 0.6, 0.4, 0.2, 0),
    c(0, 0.5, 1, 0.25, 0.75, 0.1)
  )
  raw <- normalized
  raw[, 1:4] <- sweep(
    sweep(raw[, 1:4], 2, c(10, 20, 30, 40), `*`),
    2, c(-5, 100, 1000, -50), `+`
  )

  prepared_raw <- magp:::.magp2d_prepare_fit_inputs(raw, y, q = 4)
  prepared_normalized <- magp:::.magp2d_prepare_fit_inputs(
    normalized, y, q = 4
  )
  expect_equal(prepared_raw$X, prepared_normalized$X, tolerance = 1e-14)
  expect_true(all(prepared_raw$quantity_scaling$applied))

  prepared_prediction <- magp:::.magp2d_prepare_prediction_inputs(
    raw[1:2, , drop = FALSE],
    q = 4,
    training_names = colnames(prepared_raw$X),
    quantity_scaling = prepared_raw$quantity_scaling
  )
  expect_equal(
    prepared_prediction,
    normalized[1:2, , drop = FALSE],
    tolerance = 1e-14
  )

  expect_warning(
    raw_fit <- magp2d_fit(raw, y, q = 4, seed = 1, maxeval = 1),
    "did not report convergence"
  )
  expect_warning(
    normalized_fit <- magp2d_fit(
      normalized, y, q = 4, seed = 1, maxeval = 1
    ),
    "did not report convergence"
  )
  expect_equal(raw_fit$X, normalized_fit$X, tolerance = 1e-14)
  expect_equal(raw_fit$solution, normalized_fit$solution, tolerance = 1e-12)
  expect_equal(
    predict(raw_fit, raw[1:2, , drop = FALSE]),
    predict(normalized_fit, normalized[1:2, , drop = FALSE]),
    tolerance = 1e-12
  )

  outside <- raw[1:2, , drop = FALSE]
  outside[1, 1] <- max(raw[, 1]) + 1
  expect_error(
    magp:::.magp2d_prepare_prediction_inputs(
      outside,
      q = 4,
      training_names = colnames(prepared_raw$X),
      quantity_scaling = prepared_raw$quantity_scaling
    ),
    "training range"
  )

  constant <- raw
  constant[, 1] <- 5
  prepared_constant <- magp:::.magp2d_prepare_fit_inputs(
    constant, y, q = 4
  )
  expect_true(prepared_constant$quantity_scaling$constant[1])
  expect_equal(
    unname(prepared_constant$X[, 1]),
    rep(0, nrow(constant))
  )
  changed_constant <- constant[1:2, , drop = FALSE]
  changed_constant[1, 1] <- 6
  expect_error(
    magp:::.magp2d_prepare_prediction_inputs(
      changed_constant,
      q = 4,
      training_names = colnames(prepared_constant$X),
      quantity_scaling = prepared_constant$quantity_scaling
    ),
    "constant in the training data"
  )

  wrong_halves <- X[, c(5:8, 1:4)]
  expect_error(
    magp2d_fit(wrong_halves, y),
    "permutation"
  )

  # Older fitted objects may not contain scaling metadata.
  fit <- reference_model(dat)
  bad_prediction <- dat$X_test[1:2, , drop = FALSE]
  bad_prediction[1, 2] <- 2
  expect_error(
    predict(fit, bad_prediction),
    "scaled to \\[0, 1\\]"
  )
})

test_that("prediction checks and aligns named columns", {
  dat <- example_data()
  fit <- reference_model(dat)

  newdata <- dat$test[1:4, dat$cols]
  expected <- predict(fit, newdata)
  shuffled <- newdata[, rev(dat$cols)]
  expect_equal(predict(fit, shuffled), expected)
  expected_uq <- predict(fit, newdata, se.fit = TRUE, type = "response")
  shuffled_uq <- predict(fit, shuffled, se.fit = TRUE, type = "response")
  expect_equal(shuffled_uq, expected_uq)

  wrong_names <- newdata
  names(wrong_names)[1] <- "wrong"
  expect_error(predict(fit, wrong_names), "column names")

  invalid_order <- newdata
  invalid_order[1, 5:8] <- c(1, 1, 3, 4)
  expect_error(predict(fit, invalid_order), "permutation")

  missing_value <- newdata
  missing_value[1, 1] <- NA_real_
  expect_error(predict(fit, missing_value), "finite")
  expect_error(magp:::predict.magp2d(list(), newdata), "class 'magp2d'")
})

test_that("predictive uncertainty matches an independent fixture", {
  fit <- uncertainty_fixture()
  novel <- matrix(
    c(0.25, 0.5, 0.75, 2, 3, 1),
    nrow = 1L,
    dimnames = list(NULL, colnames(fit$X))
  )

  script <- predict(fit, novel, se.fit = TRUE, type = "script")
  response <- predict(fit, novel, se.fit = TRUE, type = "response")
  latent <- predict(fit, novel, se.fit = TRUE, type = "latent")

  expect_named(script, c("fit", "se.fit", "variance", "type"))
  expect_equal(script$fit, 1.6040747686323422, tolerance = 1e-12)
  expect_equal(script$variance, 0.6120151370172487, tolerance = 1e-12)
  expect_equal(response$fit, script$fit, tolerance = 1e-14)
  expect_equal(response$variance, script$variance, tolerance = 1e-14)
  expect_equal(latent$fit, script$fit, tolerance = 1e-14)
  expect_equal(latent$variance, 0.5620151370172487, tolerance = 1e-12)
  expect_equal(response$variance - latent$variance, fit$tau, tolerance = 1e-14)

  for (out in list(script, response, latent)) {
    expect_true(all(is.finite(unlist(out[1:3]))))
    expect_true(all(out$variance >= 0))
    expect_equal(out$se.fit^2, out$variance, tolerance = 1e-14)
  }
})

test_that("training-point uncertainty distinguishes nugget conventions", {
  fit <- uncertainty_fixture()
  training_point <- fit$X[1, , drop = FALSE]

  script <- predict(fit, training_point, se.fit = TRUE, type = "script")
  response <- predict(fit, training_point, se.fit = TRUE, type = "response")
  latent <- predict(fit, training_point, se.fit = TRUE, type = "latent")

  expect_equal(script$fit, fit$y[1], tolerance = 1e-12)
  expect_equal(script$variance, 0, tolerance = 1e-10)
  expect_equal(response$fit, 1.9688319978591691, tolerance = 1e-12)
  expect_equal(latent$fit, response$fit, tolerance = 1e-14)
  expect_equal(response$variance, 0.09546958409798759, tolerance = 1e-12)
  expect_equal(latent$variance, 0.04546958409798758, tolerance = 1e-12)
  expect_equal(response$variance - latent$variance, fit$tau, tolerance = 1e-14)
})

test_that("zero nugget interpolates under every prediction convention", {
  fit <- uncertainty_fixture(tau = 0)
  training_points <- fit$X[1:2, , drop = FALSE]

  for (type in c("script", "response", "latent")) {
    out <- predict(fit, training_points, se.fit = TRUE, type = type)
    expect_equal(out$fit, fit$y[1:2], tolerance = 1e-12)
    expect_true(all(out$variance >= 0))
    expect_equal(out$variance, c(0, 0), tolerance = 1e-10)
  }
})

test_that("five-decimal script matching remains explicit", {
  fit <- uncertainty_fixture(tau = 0)
  below_boundary <- fit$X[1, , drop = FALSE]
  above_boundary <- fit$X[1, , drop = FALSE]
  below_boundary[1, 1] <- below_boundary[1, 1] + 4e-6
  above_boundary[1, 1] <- above_boundary[1, 1] + 6e-6

  below_script_fit <- predict(fit, below_boundary, type = "script")
  below_script_uq <- predict(
    fit, below_boundary, se.fit = TRUE, type = "script"
  )
  below_response <- predict(
    fit, below_boundary, se.fit = TRUE, type = "response"
  )
  above_script <- predict(
    fit, above_boundary, se.fit = TRUE, type = "script"
  )
  above_response <- predict(
    fit, above_boundary, se.fit = TRUE, type = "response"
  )

  expect_false(isTRUE(all.equal(
    below_response$fit, below_script_fit, tolerance = 1e-12
  )))
  expect_equal(below_script_uq$fit, below_script_fit, tolerance = 1e-14)
  expect_equal(below_script_uq$variance, 0, tolerance = 1e-10)
  expect_equal(above_script$fit, above_response$fit, tolerance = 1e-14)
  expect_equal(
    above_script$variance, above_response$variance, tolerance = 1e-14
  )
})

test_that("duplicate inputs do not silently return negative uncertainty", {
  base <- uncertainty_fixture()
  duplicate_fit <- base
  duplicate_fit$X <- rbind(base$X[1, ], base$X[1, ], base$X[2, ])
  duplicate_fit$y <- c(2, 2.5, -1)
  duplicate_fit$n <- 3L
  point <- duplicate_fit$X[1, , drop = FALSE]

  expect_error(
    predict(duplicate_fit, point, se.fit = TRUE, type = "script"),
    "negative predictive variance"
  )
  for (type in c("response", "latent")) {
    out <- predict(duplicate_fit, point, se.fit = TRUE, type = type)
    expect_true(all(is.finite(unlist(out[1:3]))))
    expect_true(all(out$variance >= 0))
  }
})

test_that("prediction uncertainty arguments are validated", {
  fit <- uncertainty_fixture()
  novel <- fit$X[1, , drop = FALSE]

  expect_error(predict(fit, novel, se.fit = NA), "se.fit")
  expect_error(predict(fit, novel, se.fit = 1), "se.fit")
  expect_error(predict(fit, novel, se.fit = c(TRUE, FALSE)), "se.fit")
  expect_error(predict(fit, novel, type = "unknown"), "arg")
})

test_that("uncertainty agrees across cached and reconstructed model state", {
  uncached <- uncertainty_fixture()
  state <- magp:::.magp2d_model_state(
    uncached$X, uncached$y, uncached$solution, uncached$q, uncached$tau
  )
  cached <- uncached
  cached$inverse_covariance <- state$inverse_covariance
  cached$mean <- state$mean
  cached$cholesky <- state$cholesky
  missing_cholesky <- cached
  missing_cholesky$cholesky <- NULL
  novel <- matrix(
    c(0.25, 0.5, 0.75, 2, 3, 1),
    nrow = 1L,
    dimnames = list(NULL, colnames(uncached$X))
  )

  for (type in c("response", "latent")) {
    expected <- predict(uncached, novel, se.fit = TRUE, type = type)
    expect_equal(
      predict(cached, novel, se.fit = TRUE, type = type), expected
    )
    expect_equal(
      predict(missing_cholesky, novel, se.fit = TRUE, type = type),
      expected
    )
  }
})

test_that("optimizer limits are surfaced", {
  dat <- example_data()
  expect_warning(
    fit <- magp2d_fit(
      dat$X[1:6, , drop = FALSE], dat$y[1:6],
      q = 4, seed = 1, maxeval = 1
    ),
    "did not report convergence"
  )
  expect_false(fit$converged)
})

test_that("RMSE rejects malformed values", {
  expect_equal(magp2d_rmse(c(1, 2), c(1, 4)), sqrt(2))
  expect_error(magp2d_rmse(1:2, 1:3), "same length")
  expect_error(magp2d_rmse(numeric(), numeric()), "must not be empty")
  expect_error(magp2d_rmse(c(1, NA), c(1, 2)), "finite")
  expect_error(magp2d_rmse(factor(c("a", "b")), 1:2), "numeric")
})
