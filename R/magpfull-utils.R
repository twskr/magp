# Internal helpers for the full-mapping MaGP model.

#' Number of fitted covariance parameters in a full-MaGP model
#'
#' @noRd
.magpfull_npar <- function(q) {
  as.integer(2L * q + q * (q - 1L) / 2L)
}

#' Map an order level to its full-mapping latent vector
#'
#' @noRd
.magpfull_omap <- function(l, delta, q) {
  mapping_dimension <- q - 1L
  expected <- q * (q - 1L) / 2L
  if (length(delta) != expected) {
    stop(
      "full-mapping delta must have length q*(q-1)/2",
      call. = FALSE
    )
  }
  triangle <- matrix(0, mapping_dimension, mapping_dimension)
  upper_positions <- which(upper.tri(triangle, diag = TRUE), arr.ind = TRUE)
  triangle[upper_positions] <- delta
  mapping <- rbind(rep(0, mapping_dimension), t(triangle))
  mapping[l, , drop = FALSE]
}

#' Locate a full-mapping parameter in the triangular mapping
#'
#' @noRd
.magpfull_delta_position <- function(s, q) {
  for (level in 2L:q) {
    lower <- (level - 1L) * (level - 2L) / 2L
    upper <- level * (level - 1L) / 2L
    if (lower < s && s <= upper) {
      return(c(level = level, coordinate = s - lower))
    }
  }
  stop("invalid full-mapping parameter index", call. = FALSE)
}
