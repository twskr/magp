full_example_data <- function() {
  train <- read.table(
    system.file("extdata", "example_train.txt", package = "magp"),
    header = TRUE
  )
  test <- read.table(
    system.file("extdata", "example_test.txt", package = "magp"),
    header = TRUE
  )
  columns <- c("A", "B", "C", "D", "a", "b", "c", "d")
  list(
    train = train,
    test = test,
    columns = columns,
    X = as.matrix(train[, columns]),
    y = train$y,
    X_test = as.matrix(test[, columns]),
    y_test = test$y
  )
}

full_reference_solution <- c(
  272.800857744284, 378.335163658610, 577.059319563649,
  909.066070791303, 202.069716717763, 898.439430888460,
  944.700330812723, 660.957723779057, -0.00180012951560266,
  -0.00140792885480384, 0.00213475584125560, -0.000444960944009809,
  -0.000680506855395875, -0.000932213757011885
)

full_reference_model <- function(data, tau = 0.001) {
  structure(
    list(
      mapping = "full",
      mapping_dimension = 3L,
      q = 4L,
      n = nrow(data$X),
      tau = tau,
      npar = 14L,
      X = data$X,
      y = data$y,
      solution = full_reference_solution
    ),
    class = "magpfull"
  )
}

test_that("full mapping follows the triangular layout", {
  expect_equal(
    magp:::.magpfull_omap(1:4, 1:6, q = 4),
    rbind(
      c(0, 0, 0),
      c(1, 0, 0),
      c(2, 3, 0),
      c(4, 5, 6)
    )
  )
  expect_equal(
    vapply(
      1:6,
      function(s) magp:::.magpfull_delta_position(s, 4)["level"],
      numeric(1)
    ),
    c(2, 3, 3, 4, 4, 4)
  )
  expect_equal(
    vapply(
      1:6,
      function(s) {
        magp:::.magpfull_delta_position(s, 4)["coordinate"]
      },
      numeric(1)
    ),
    c(1, 1, 2, 1, 2, 3)
  )
})

test_that("full-mapping analytical gradients match finite differences", {
  data <- full_example_data()
  X <- data$X[1:8, , drop = FALSE]
  y <- data$y[1:8]
  evaluate <- magp:::.magpfull_objective_factory(X, y, 4, 0.001)
  parameter <- c(
    120, 90, 110, 80,
    4, 7, 3, 5,
    0.2, -0.3, 0.4, -0.1, 0.25, -0.2
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

  expect_length(analytical, 14L)
  expect_lt(max(scaled_error), 1e-6)
})

test_that("full fit reproduces the stored numerical reference", {
  data <- full_example_data()
  fit <- magpfull_fit(
    data$X, data$y,
    q = 4, tau = 0.001, seed = 1, maxeval = 500
  )

  expect_s3_class(fit, "magpfull")
  expect_true(fit$converged)
  expect_identical(fit$mapping, "full")
  expect_equal(fit$mapping_dimension, 3L)
  expect_equal(fit$npar, 14L)
  expect_equal(fit$solution, full_reference_solution, tolerance = 1e-7)
  expect_equal(fit$objective, 237.56338710779033, tolerance = 1e-8)
  expect_named(fit$parameters, c("sigma2", "theta", "delta"))
  expect_length(fit$parameters$sigma2, 4L)
  expect_length(fit$parameters$theta, 4L)
  expect_length(fit$parameters$delta, 6L)
  expect_length(fit$initial, 14L)

  prediction <- predict(fit, data$X_test)
  expect_identical(
    prediction,
    predict(fit, data$X_test, se.fit = FALSE, type = "script")
  )
  expect_equal(
    magp2d_rmse(prediction, data$y_test),
    11.388894042229163,
    tolerance = 1e-10
  )
})

test_that("full fit accepts a y column and infers q", {
  data <- full_example_data()
  expect_warning(
    fit <- magpfull_fit(data$train, seed = 1, maxeval = 1),
    "did not report convergence"
  )

  expect_equal(fit$q, 4L)
  expect_equal(fit$y, data$train$y)
  expect_equal(colnames(fit$X), data$columns)
  expect_equal(fit$input_roles$quantity, c("A", "B", "C", "D"))
  expect_equal(fit$input_roles$sequence, c("a", "b", "c", "d"))
  expect_named(fit$input_roles, c("quantity", "sequence"))
  expect_false(any(fit$quantity_scaling$applied))
  expect_equal(
    predict(fit, data$test[1:3, ]),
    predict(fit, data$X_test[1:3, , drop = FALSE])
  )
})

test_that("full prediction implements response and latent uncertainty", {
  data <- full_example_data()
  fit <- full_reference_model(data)
  novel <- data$X_test[1:3, , drop = FALSE]

  script <- predict(fit, novel, se.fit = TRUE, type = "script")
  response <- predict(fit, novel, se.fit = TRUE, type = "response")
  latent <- predict(fit, novel, se.fit = TRUE, type = "latent")

  expect_named(response, c("fit", "se.fit", "variance", "type"))
  expect_equal(script$fit, response$fit, tolerance = 1e-14)
  expect_equal(script$variance, response$variance, tolerance = 1e-14)
  expect_equal(response$fit, latent$fit, tolerance = 1e-14)
  expect_equal(
    response$variance - latent$variance,
    rep(fit$tau, nrow(novel)),
    tolerance = 1e-9
  )
  expect_true(all(response$variance >= 0))
  expect_true(all(latent$variance >= 0))
  expect_equal(response$se.fit^2, response$variance, tolerance = 1e-14)
  expect_equal(latent$se.fit^2, latent$variance, tolerance = 1e-14)
})

test_that("full script mode interpolates a unique exact training row", {
  data <- full_example_data()
  fit <- full_reference_model(data)
  training_row <- data$X[1, , drop = FALSE]

  script <- predict(
    fit, training_row, se.fit = TRUE, type = "script"
  )
  response <- predict(
    fit, training_row, se.fit = TRUE, type = "response"
  )
  latent <- predict(
    fit, training_row, se.fit = TRUE, type = "latent"
  )

  expect_equal(script$fit, data$y[1], tolerance = 1e-10)
  expect_equal(script$variance, 0, tolerance = 1e-8)
  expect_equal(response$fit, latent$fit, tolerance = 1e-14)
  expect_equal(
    response$variance - latent$variance,
    fit$tau,
    tolerance = 1e-9
  )
})

test_that("full prediction aligns columns and reconstructs cached state", {
  data <- full_example_data()
  uncached <- full_reference_model(data)
  state <- magp:::.magpfull_model_state(
    uncached$X,
    uncached$y,
    uncached$solution,
    uncached$q,
    uncached$tau
  )
  cached <- uncached
  cached$mean <- state$mean
  cached$cholesky <- state$cholesky
  cached$inverse_covariance <- state$inverse_covariance
  newdata <- data$test[1:4, data$columns]
  shuffled <- newdata[, rev(data$columns)]

  expected <- predict(
    uncached, newdata, se.fit = TRUE, type = "response"
  )
  expect_equal(
    predict(cached, newdata, se.fit = TRUE, type = "response"),
    expected
  )
  expect_equal(
    predict(cached, shuffled, se.fit = TRUE, type = "response"),
    expected
  )
})

test_that("full-mapping inputs and prediction arguments are validated", {
  data <- full_example_data()
  X <- data$X[1:6, , drop = FALSE]
  y <- data$y[1:6]
  bad_order <- X
  bad_order[1, 5:8] <- c(1, 1, 3, 4)

  expect_error(magpfull_fit(X, y, q = 1), "q")
  expect_error(magpfull_fit(bad_order, y, q = 4), "permutation")
  expect_error(magpfull_fit(X, y, q = 4, tau = -1), "tau")
  expect_error(
    magpfull_fit(X, y, q = 4, lb_delta = 1, ub_delta = -1),
    "lower parameter bound"
  )

  fit <- full_reference_model(data)
  point <- data$X_test[1, , drop = FALSE]
  expect_error(predict(fit, point, se.fit = NA), "se.fit")
  expect_error(predict(fit, point, type = "unknown"), "arg")
  expect_error(
    magp:::predict.magpfull(list(), point),
    "class 'magpfull'"
  )
})
