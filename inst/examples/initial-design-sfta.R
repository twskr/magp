# Run after installing magp >= 0.12.0.
library(magp)
stopifnot(packageVersion("magp") >= "0.12.0")

# Compare random and SFTA sequence designs using the same quantitative design.
random_design <- magp_initial_design(
  n = 16, q = 4, sequence_method = "random", seed = 1,
  quantity_maxit = 1000, alignment_maxit = 1000
)
sfta_design <- magp_initial_design(
  n = 16, q = 4, sequence_method = "sfta", seed = 1,
  sequence_maxit = 2000, quantity_maxit = 1000, alignment_maxit = 1000,
  sfta_control = list(nstarts = 5, ncalibrate = 200, nrounds = 10)
)
print(rbind(random = random_design$criteria, sfta = sfta_design$criteria))
stopifnot(identical(random_design$quantity_search$quantity,
                    sfta_design$quantity_search$quantity))
stopifnot(sfta_design$sequence_search$criterion <=
            sfta_design$sequence_search$initial_criterion)
print(sfta_design$design)
print(sfta_design$sequence_search$sfta$thresholds)

# The first q columns contain quantitative inputs. The last q columns give
# component positions, with each row a permutation of 1:q.

# Use SFTA for the initial design in a Bayesian optimization run.
objective <- function(quantity_1, quantity_2, quantity_3,
                      sequence_1, sequence_2, sequence_3) {
  x <- c(quantity_1, quantity_2, quantity_3)
  o <- c(sequence_1, sequence_2, sequence_3)
  sum((x - c(0.2, 0.6, 0.8))^2) + 0.01 * sum((o - c(1, 3, 2))^2)
}
bo <- magp_bayes_optimize_from_scratch(
  FUN = objective, n_initial = 8, q = 3, n_iter = 1,
  direction = "minimize", model = "2d", seed = 4,
  design_control = list(
    sequence_method = "sfta", sequence_maxit = 500,
    quantity_maxit = 500, alignment_maxit = 500,
    sfta_control = list(nstarts = 3, ncalibrate = 100)
  ),
  fit_control = list(maxeval = 500),
  acquisition_control = list(max_sequences = 6, n_starts = 2, maxit = 30),
  verbose = FALSE
)
print(bo$best_point)
print(bo$best_value)
print(bo$history)
stopifnot(bo$initial_design$sequence_method == "sfta",
          bo$iterations_completed == 1L)

# SFTA changes the initial sequence design, not the later EI search.
