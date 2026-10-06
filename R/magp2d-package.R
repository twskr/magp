#' magp: Mapping-Based Additive Gaussian Process Models
#'
#' Fits additive Gaussian process models for experiments in which each
#' component has both a quantitative level and a position in a sequence.
#'
#' For `q` components, the first `q` input columns contain quantitative levels
#' and the next `q` columns contain sequence positions. Every sequence row must
#' be a permutation of `1:q`. Quantitative columns outside `[0, 1]` are scaled
#' during fitting, and the same transformation is used for new data.
#'
#' [magp2d_fit()] represents the sequence positions in two latent dimensions.
#' [magpfull_fit()] uses `q - 1` latent dimensions. Both fitted model classes
#' support [stats::predict()] for point predictions and plug-in predictive
#' uncertainty. The computational kernels for covariance matrices, analytical
#' gradients, and cross-covariances are implemented in C++ with `Rcpp`.
#'
#' [magp_initial_design()] constructs a space-filling quantitative-sequence
#' design before responses are collected. It combines a Latin hypercube with
#' sequence permutations generated randomly or improved with simulated
#' annealing or space-filling threshold accepting. Joint alignment preserves
#' both component designs.
#'
#' [magp_expected_improvement()] evaluates improvement using latent predictive
#' uncertainty. [magp_next_point()] searches quantitative bounds and sequence
#' permutations for the next experiment, while [magp_bayes_optimize()] runs
#' the sequential fitting and evaluation loop from completed experiments.
#' [magp_bayes_optimize_from_scratch()] first generates and evaluates an
#' initial design, then continues through the same sequential loop.
#'
#' @useDynLib magp, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @keywords internal
"_PACKAGE"
