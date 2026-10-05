# SFTA for the complete sequence-design criterion.

.magp_sfta_control <- function(control, n, maxit) {
  if (!is.list(control) || (length(control) &&
      (is.null(names(control)) || anyNA(names(control)) ||
       any(!nzchar(names(control))) || anyDuplicated(names(control))))) {
    stop("sfta_control must be a uniquely named list", call. = FALSE)
  }
  defaults <- list(
    nstarts = 5L, ncalibrate = 200L, nrounds = min(10L, maxit),
    max_proposals = min(.Machine$integer.max, max(1000, 50 * as.double(n)))
  )
  unknown <- setdiff(names(control), names(defaults))
  if (length(unknown)) {
    stop("unknown sfta_control argument(s): ", paste(unknown, collapse = ", "),
         call. = FALSE)
  }
  defaults[names(control)] <- control
  for (key in names(defaults)) {
    .magp2d_validate_scalar(
      defaults[[key]], paste0("sfta_control$", key), lower = 1,
      upper = .Machine$integer.max, integer = TRUE
    )
    defaults[[key]] <- as.integer(defaults[[key]])
  }
  if (defaults$nrounds > maxit) {
    stop("sfta_control$nrounds must not exceed maxit (sequence_maxit)",
         call. = FALSE)
  }
  defaults
}

# Bounded Phase-I sampler: a row is accepted with probability H_min / q.
# Once the finite space is exhausted, or the budget is reached, fill remaining
# rows uniformly. This explicitly allows repeated sequences if n > q!.
.magp_sfta_space_filling_orders <- function(n, q, max_proposals) {
  orders <- matrix(NA_integer_, n, q)
  orders[1L, ] <- sample.int(q)
  filled <- 1L
  attempts <- 1L
  finite_space <- if (lgamma(q + 1) <= log(n) + 1e-12) {
    round(exp(lgamma(q + 1)))
  } else {
    Inf
  }
  while (filled < n && attempts < max_proposals && filled < finite_space) {
    candidate <- sample.int(q)
    attempts <- attempts + 1L
    hmin <- min(rowSums(sweep(
      orders[seq_len(filled), , drop = FALSE], 2L, candidate, FUN = "!="
    )))
    if (stats::runif(1L) < hmin / q) {
      filled <- filled + 1L
      orders[filled, ] <- candidate
    }
  }
  fallback <- n - filled
  if (fallback > 0L) {
    for (i in seq.int(filled + 1L, n)) orders[i, ] <- sample.int(q)
  }
  list(orders = orders, attempts = attempts, fallbacks = fallback)
}

.magp_sfta_sequence_search <- function(initial_orders, weights, p, maxit,
                                      control) {
  n <- nrow(initial_orders)
  q <- ncol(initial_orders)
  best_orders <- initial_orders
  initial_criterion <- .magp_sequence_criterion_orders(initial_orders, weights, p)
  best_criterion <- initial_criterion
  attempts <- fallbacks <- integer(control$nstarts)
  phase1_values <- numeric(control$nstarts)
  for (i in seq_len(control$nstarts)) {
    sampled <- .magp_sfta_space_filling_orders(n, q, control$max_proposals)
    value <- .magp_sequence_criterion_orders(sampled$orders, weights, p)
    attempts[i] <- sampled$attempts
    fallbacks[i] <- sampled$fallbacks
    phase1_values[i] <- value
    if (value < best_criterion) {
      best_orders <- sampled$orders
      best_criterion <- value
    }
  }
  phase1_criterion <- best_criterion
  fit <- magp_sfta_search_cpp(
    best_orders, weights[["pair_balance"]], weights[["space_filling"]],
    p, maxit, control$ncalibrate, control$nrounds
  )
  # Independently recompute the public criterion rather than expose cached sums.
  final_criterion <- .magp_sequence_criterion_orders(fit$orders, weights, p)
  if (final_criterion > best_criterion) {
    fit$orders <- best_orders
    final_criterion <- best_criterion
  }
  list(
    orders = fit$orders, criterion = final_criterion,
    initial_criterion = initial_criterion,
    diagnostics = list(
      controls = control, phase1_criterion = phase1_criterion,
      phase1_values = phase1_values, phase1_attempts = attempts,
      phase1_fallbacks = fallbacks, thresholds = fit$thresholds,
      accepted = fit$accepted, best_history = fit$best_history,
      iterations = maxit,
      evaluations = fit$evaluations + control$nstarts + 2
    )
  )
}
