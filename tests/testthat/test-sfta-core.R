test_that("incremental SFTA scores agree with the full public criterion", {
  cases <- expand.grid(
    n = c(2L, 6L, 8L, 40L), q = c(3L, 4L, 8L),
    p = c(1L, 15L), pair_weight = c(0, 0.2, 1)
  )
  for (i in seq_len(nrow(cases))) {
    settings <- cases[i, ]
    set.seed(7000L + i)
    orders <- t(replicate(settings$n, sample.int(settings$q)))
    fit <- magp:::magp_sfta_search_cpp(
      orders, settings$pair_weight, 1 - settings$pair_weight,
      settings$p, 250L, 31L, 7L
    )
    to_positions <- function(x) t(apply(x, 1L, order))
    expect_equal(
      fit$initial_criterion,
      magp_sequence_criterion(
        to_positions(orders), pair_weight = settings$pair_weight,
        space_weight = 1 - settings$pair_weight, p = settings$p
      ), tolerance = 1e-12
    )
    expect_equal(
      fit$criterion,
      magp_sequence_criterion(
        to_positions(fit$orders), pair_weight = settings$pair_weight,
        space_weight = 1 - settings$pair_weight, p = settings$p
      ), tolerance = 1e-12
    )
    expect_lte(fit$criterion, fit$initial_criterion)
    expect_true(all(apply(fit$orders, 1L, function(row) {
      identical(sort(row), seq_len(settings$q))
    })))
    expect_equal(fit$iterations, 250L)
    expect_equal(fit$evaluations, 282)
    expect_length(fit$best_history, 8L)
    expect_equal(fit$best_history[1L], fit$initial_criterion)
    expect_equal(tail(fit$best_history, 1L), fit$criterion)
    expect_true(all(diff(fit$best_history) <= 0))
    expect_true(all(diff(fit$thresholds) <= 0))
    expect_equal(tail(fit$thresholds, 1L), 0)
    expect_equal(sum(fit$accepted_per_round), fit$accepted)
  }
})

test_that("long SFTA runs do not accumulate incremental criterion drift", {
  for (settings in list(c(n = 8L, q = 3L), c(n = 40L, q = 8L))) {
    set.seed(876)
    orders <- t(replicate(settings[["n"]], sample.int(settings[["q"]])))
    fit <- magp:::magp_sfta_search_cpp(
      orders, 0.2, 0.8, 15L, 100000L, 1000L, 20L
    )
    positions <- t(apply(fit$orders, 1L, order))
    expect_equal(fit$criterion, magp_sequence_criterion(positions),
                 tolerance = 1e-12)
    expect_lte(fit$criterion, fit$initial_criterion)
    expect_true(all(diff(fit$best_history) <= 0))
    expect_equal(fit$iterations, 100000L)
    expect_equal(fit$evaluations, 101001)
  }
})

test_that("single-round and flat-objective SFTA searches are bounded", {
  orders <- rbind(1:3, 1:3)
  set.seed(916)
  fit <- magp:::magp_sfta_search_cpp(orders, 0, 1, 15L, 17L, 1L, 1L)
  expect_identical(fit$thresholds, 0)
  expect_length(fit$best_history, 2L)
  expect_equal(fit$iterations, 17L)
  expect_equal(fit$evaluations, 19)
  expect_lte(fit$criterion, fit$initial_criterion)

  # Equal current/proposed scores must not pass the final strict greedy test.
  separated <- rbind(1:3, 3:1)
  set.seed(111)
  fit <- magp:::magp_sfta_search_cpp(separated, 0, 1, 15L, 100L, 20L, 1L)
  expect_equal(fit$accepted, 1L)
  expect_equal(fit$criterion, (0.25^15)^(1 / 15))
})

test_that("the SFTA core is reproducible and never mutates its input", {
  orders <- rbind(1:4, 4:1, c(2L, 4L, 1L, 3L))
  original <- orders
  set.seed(532)
  first <- magp:::magp_sfta_search_cpp(orders, .2, .8, 15L, 101L, 23L, 7L)
  state_after <- .Random.seed
  set.seed(532)
  second <- magp:::magp_sfta_search_cpp(orders, .2, .8, 15L, 101L, 23L, 7L)
  expect_identical(first, second)
  expect_identical(orders, original)
  expect_identical(.Random.seed, state_after)
  set.seed(532)
  invisible(runif(3L * (101L + 23L)))
  expect_identical(.Random.seed, state_after)
})

test_that("the internal SFTA core validates its dimensions and permutations", {
  valid <- rbind(1:3, 3:1)
  search <- function(orders = valid, pair_weight = .2, space_weight = .8,
                     p = 15L, maxit = 10L, ncalibrate = 2L, nrounds = 2L) {
    magp:::magp_sfta_search_cpp(
      orders, pair_weight, space_weight, p, maxit, ncalibrate, nrounds
    )
  }
  expect_error(search(matrix(1:3, nrow = 1L)), "two rows")
  expect_error(search(rbind(1:2, 2:1)), "three columns")
  for (replacement in c(0L, 4L, NA_integer_, 2L)) {
    invalid <- valid
    invalid[1L, 1L] <- replacement
    expect_error(search(invalid), "permutation")
  }
  expect_error(search(pair_weight = -1), "weights")
  expect_error(search(space_weight = Inf), "weights")
  expect_error(search(pair_weight = NaN), "weights")
  expect_error(search(pair_weight = 0, space_weight = 0), "positive sum")
  expect_error(search(p = 0L), "positive")
  expect_error(search(maxit = 0L), "positive")
  expect_error(search(ncalibrate = 0L), "positive")
  expect_error(search(nrounds = 0L), "positive")
  expect_error(search(nrounds = 11L), "must not exceed")
  expect_error(search(p = NA_integer_), "positive")
})
