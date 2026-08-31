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
#' @useDynLib magp, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @keywords internal
"_PACKAGE"
