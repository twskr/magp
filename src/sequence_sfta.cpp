#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>

namespace {

// The sequence criterion is evaluated on execution-order rows.
// Integer histograms are updated exactly. Rebuilding the two positive sums
// from these histograms avoids cancellation and accumulated score drift.
class SequenceCriterionState {
 public:
  const int n;
  const int q;
  std::vector<int> orders;

  SequenceCriterionState(const Rcpp::IntegerMatrix& input,
                         double pair_weight, double space_weight, int p)
      : n(input.nrow()), q(input.ncol()),
        orders(static_cast<std::size_t>(n) * q),
        pair_counts(static_cast<std::size_t>(q) * q, 0),
        hamming(static_cast<std::size_t>(n) * (n - 1) / 2, 0),
        pair_histogram(static_cast<std::size_t>(n) + 1, 0.0),
        hamming_histogram(static_cast<std::size_t>(q) + 1, 0.0),
        pair_power(static_cast<std::size_t>(n) + 1),
        hamming_power(static_cast<std::size_t>(q) + 1),
        pair_weight_(pair_weight), space_weight_(space_weight),
        inverse_p(1.0 / p) {
    for (int i = 0; i < n; ++i) {
      if (i % 128 == 0) Rcpp::checkUserInterrupt();
      for (int j = 0; j < q; ++j) {
        orders[row_index(i, j)] = input(i, j) - 1;
      }
      for (int j = 0; j < q - 1; ++j) {
        ++pair_counts[edge_index(orders[row_index(i, j)],
                                 orders[row_index(i, j + 1)])];
      }
      for (int other = 0; other < i; ++other) {
        int distance = 0;
        for (int j = 0; j < q; ++j) {
          distance += orders[row_index(i, j)] != orders[row_index(other, j)];
        }
        hamming[distance_index(i, other)] = distance;
        ++hamming_histogram[distance];
      }
    }
    for (int a = 0; a < q; ++a) {
      for (int b = 0; b < q; ++b) {
        if (a != b) ++pair_histogram[pair_counts[edge_index(a, b)]];
      }
    }
    for (int count = 0; count <= n; ++count) {
      pair_power[count] = std::pow(static_cast<double>(count) + 1.0, -p);
    }
    for (int distance = 0; distance <= q; ++distance) {
      hamming_power[distance] =
          std::pow(static_cast<double>(distance) + 1.0, -p);
    }
  }

  double score() const {
    long double pair_sum = 0.0L;
    long double hamming_sum = 0.0L;
    for (int count = 0; count <= n; ++count) {
      pair_sum += static_cast<long double>(pair_histogram[count]) *
          pair_power[count];
    }
    for (int distance = 0; distance <= q; ++distance) {
      hamming_sum += static_cast<long double>(hamming_histogram[distance]) *
          hamming_power[distance];
    }
    const double combined = static_cast<double>(
        pair_weight_ * pair_sum + space_weight_ * hamming_sum);
    return std::pow(combined, inverse_p);
  }

  // This operation is an involution: calling it again with the same indices
  // restores every integer count/distance and the original execution order.
  void swap_columns(int row, int a, int b) {
    int edges[4] = {a - 1, a, b - 1, b};
    std::sort(edges, edges + 4);
    const int n_edges = static_cast<int>(std::unique(edges, edges + 4) - edges);
    for (int k = 0; k < n_edges; ++k) {
      const int edge = edges[k];
      if (edge >= 0 && edge < q - 1) {
        change_pair(orders[row_index(row, edge)],
                    orders[row_index(row, edge + 1)], -1);
      }
    }

    const int old_a = orders[row_index(row, a)];
    const int old_b = orders[row_index(row, b)];
    for (int other = 0; other < n; ++other) {
      if (other == row) continue;
      const int other_a = orders[row_index(other, a)];
      const int other_b = orders[row_index(other, b)];
      const int delta = (old_b != other_a) + (old_a != other_b) -
          (old_a != other_a) - (old_b != other_b);
      if (delta != 0) {
        const std::size_t index = distance_index(row, other);
        const int old_distance = hamming[index];
        --hamming_histogram[old_distance];
        hamming[index] += delta;
        ++hamming_histogram[hamming[index]];
      }
    }
    std::swap(orders[row_index(row, a)], orders[row_index(row, b)]);

    for (int k = 0; k < n_edges; ++k) {
      const int edge = edges[k];
      if (edge >= 0 && edge < q - 1) {
        change_pair(orders[row_index(row, edge)],
                    orders[row_index(row, edge + 1)], 1);
      }
    }
  }

 private:
  std::vector<int> pair_counts;
  std::vector<int> hamming;
  std::vector<double> pair_histogram;
  std::vector<double> hamming_histogram;
  std::vector<double> pair_power;
  std::vector<double> hamming_power;
  double pair_weight_;
  double space_weight_;
  double inverse_p;

  std::size_t row_index(int row, int col) const {
    return static_cast<std::size_t>(row) * q + col;
  }

  std::size_t edge_index(int from, int to) const {
    return static_cast<std::size_t>(from) * q + to;
  }

  std::size_t distance_index(int a, int b) const {
    const int larger = std::max(a, b);
    const int smaller = std::min(a, b);
    return static_cast<std::size_t>(larger) * (larger - 1) / 2 + smaller;
  }

  void change_pair(int from, int to, int change) {
    const std::size_t index = edge_index(from, to);
    --pair_histogram[pair_counts[index]];
    pair_counts[index] += change;
    ++pair_histogram[pair_counts[index]];
  }
};

struct SequenceSwap {
  int row;
  int a;
  int b;
};

SequenceSwap random_swap(int n, int q) {
  SequenceSwap result;
  result.row = static_cast<int>(R::runif(0.0, static_cast<double>(n)));
  result.a = static_cast<int>(R::runif(0.0, static_cast<double>(q)));
  result.b = static_cast<int>(R::runif(0.0, static_cast<double>(q - 1)));
  if (result.b >= result.a) ++result.b;
  return result;
}

double empirical_quantile(const std::vector<double>& sorted, double probability) {
  const double location = (sorted.size() - 1) * probability;
  const std::size_t lower = static_cast<std::size_t>(std::floor(location));
  const std::size_t upper = static_cast<std::size_t>(std::ceil(location));
  const double fraction = location - lower;
  // R's default quantile(type = 7), including ncalibrate == 1.
  return (1.0 - fraction) * sorted[lower] + fraction * sorted[upper];
}

}  // namespace

// Threshold-acceptance phase of SFTA for the sequence initial-design criterion.
// The R wrapper supplies the space-filling Phase-I design and local RNG scope.
// [[Rcpp::export]]
Rcpp::List magp_sfta_search_cpp(Rcpp::IntegerMatrix orders, double pair_weight,
                               double space_weight, int p, int maxit,
                               int ncalibrate, int nrounds) {
  const int n = orders.nrow();
  const int q = orders.ncol();
  if (n < 2 || q < 3) Rcpp::stop("orders must have at least two rows and three columns");
  if (!std::isfinite(pair_weight) || !std::isfinite(space_weight) ||
      pair_weight < 0.0 || space_weight < 0.0 ||
      !std::isfinite(pair_weight + space_weight) ||
      pair_weight + space_weight <= 0.0) {
    Rcpp::stop("sequence weights must be finite, nonnegative, and have a positive sum");
  }
  if (p == NA_INTEGER || p < 1 || maxit == NA_INTEGER || maxit < 1 ||
      ncalibrate == NA_INTEGER || ncalibrate < 1 ||
      nrounds == NA_INTEGER || nrounds < 1 || nrounds > maxit) {
    Rcpp::stop("p, maxit, ncalibrate, and nrounds must be positive; nrounds must not exceed maxit");
  }
  std::vector<int> seen(q, -1);
  for (int i = 0; i < n; ++i) {
    if (i % 128 == 0) Rcpp::checkUserInterrupt();
    for (int j = 0; j < q; ++j) {
      const int value = orders(i, j);
      if (value == NA_INTEGER || value < 1 || value > q || seen[value - 1] == i) {
        Rcpp::stop("every row of orders must be a permutation of 1:q");
      }
      seen[value - 1] = i;
    }
  }

  SequenceCriterionState state(orders, pair_weight, space_weight, p);
  const double initial_score = state.score();
  double current_score = initial_score;
  double best_score = initial_score;
  std::vector<int> best_orders = state.orders;
  std::vector<double> calibration(ncalibrate);
  for (int i = 0; i < ncalibrate; ++i) {
    if (i % 1024 == 0) Rcpp::checkUserInterrupt();
    const SequenceSwap move = random_swap(n, q);
    state.swap_columns(move.row, move.a, move.b);
    calibration[i] = std::abs(state.score() - initial_score);
    state.swap_columns(move.row, move.a, move.b);
  }
  std::sort(calibration.begin(), calibration.end());

  Rcpp::NumericVector thresholds(nrounds);
  Rcpp::NumericVector best_history(nrounds + 1);
  Rcpp::IntegerVector accepted_per_round(nrounds);
  best_history[0] = best_score;
  int accepted = 0;
  int iterations = 0;
  for (int round = 0; round < nrounds; ++round) {
    const double probability = 0.5 * (1.0 - (round + 1.0) / nrounds);
    thresholds[round] = empirical_quantile(calibration, probability);
    // The smallest calibrated score change may be positive. Set the final
    // threshold to zero so this round accepts only strict improvements.
    if (round == nrounds - 1) thresholds[round] = 0.0;
    const int round_iterations = maxit / nrounds + (round < maxit % nrounds);
    for (int step = 0; step < round_iterations; ++step) {
      if (iterations % 1024 == 0) Rcpp::checkUserInterrupt();
      const SequenceSwap move = random_swap(n, q);
      state.swap_columns(move.row, move.a, move.b);
      const double proposal_score = state.score();
      if (proposal_score - current_score < thresholds[round]) {
        current_score = proposal_score;
        ++accepted;
        ++accepted_per_round[round];
        if (proposal_score < best_score) {
          best_score = proposal_score;
          best_orders = state.orders;
        }
      } else {
        state.swap_columns(move.row, move.a, move.b);
      }
      ++iterations;
    }
    best_history[round + 1] = best_score;
  }

  Rcpp::IntegerMatrix result(n, q);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < q; ++j) {
      result(i, j) = best_orders[static_cast<std::size_t>(i) * q + j] + 1;
    }
  }
  return Rcpp::List::create(
      Rcpp::Named("orders") = result,
      Rcpp::Named("criterion") = best_score,
      Rcpp::Named("initial_criterion") = initial_score,
      Rcpp::Named("thresholds") = thresholds,
      Rcpp::Named("accepted") = accepted,
      Rcpp::Named("accepted_per_round") = accepted_per_round,
      Rcpp::Named("evaluations") = static_cast<double>(ncalibrate) + maxit + 1.0,
      Rcpp::Named("iterations") = iterations,
      Rcpp::Named("best_history") = best_history,
      Rcpp::Named("convergence") = 0);
}
