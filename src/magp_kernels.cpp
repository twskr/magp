#include <Rcpp.h>
#include <cmath>
#include <vector>

using namespace Rcpp;

namespace {

struct MappingInfo {
  int dimension;
  std::vector<double> values;
  std::vector<int> parameter;
};

inline double mapping_value(const MappingInfo& mapping, int level, int coordinate) {
  return mapping.values[level * mapping.dimension + coordinate];
}

inline int mapping_parameter(const MappingInfo& mapping, int level, int coordinate) {
  return mapping.parameter[level * mapping.dimension + coordinate];
}

MappingInfo build_mapping(const NumericVector& parameters, int q, int mapping_type) {
  const int delta_offset = 2 * q;
  MappingInfo mapping;
  mapping.dimension = mapping_type == 0 ? 2 : q - 1;
  mapping.values.assign(q * mapping.dimension, 0.0);
  mapping.parameter.assign(q * mapping.dimension, -1);

  int s = 0;
  if (mapping_type == 0) {
    // Fix the first two positions for identifiability, then assign each
    // remaining position a pair of free coordinates.
    mapping.values[1 * mapping.dimension + 1] = parameters[delta_offset];
    mapping.parameter[1 * mapping.dimension + 1] = 0;
    s = 1;
    for (int level = 2; level < q; ++level) {
      for (int coordinate = 0; coordinate < 2; ++coordinate) {
        mapping.values[level * mapping.dimension + coordinate] =
          parameters[delta_offset + s];
        mapping.parameter[level * mapping.dimension + coordinate] = s;
        ++s;
      }
    }
  } else {
    // In the full mapping, each zero-based sequence level has the same number
    // of free coordinates as its level index.
    for (int level = 1; level < q; ++level) {
      for (int coordinate = 0; coordinate < level; ++coordinate) {
        mapping.values[level * mapping.dimension + coordinate] =
          parameters[delta_offset + s];
        mapping.parameter[level * mapping.dimension + coordinate] = s;
        ++s;
      }
    }
  }
  return mapping;
}

inline double component_base(
    const NumericMatrix& x,
    int i,
    int j,
    int h,
    int q,
    const NumericVector& parameters,
    const MappingInfo& mapping) {
  const double difference = x(i, h) - x(j, h);
  const int first_order = static_cast<int>(x(i, q + h)) - 1;
  const int second_order = static_cast<int>(x(j, q + h)) - 1;
  double order_distance = 0.0;
  for (int coordinate = 0; coordinate < mapping.dimension; ++coordinate) {
    const double delta =
      mapping_value(mapping, first_order, coordinate) -
      mapping_value(mapping, second_order, coordinate);
    order_distance += delta * delta;
  }
  return std::exp(
    -parameters[q + h] * difference * difference - order_distance
  );
}

inline double cross_component_base(
    const NumericMatrix& training,
    const NumericMatrix& prediction,
    int train_row,
    int prediction_row,
    int h,
    int q,
    const NumericVector& parameters,
    const MappingInfo& mapping) {
  const double difference =
    prediction(prediction_row, h) - training(train_row, h);
  const int prediction_order =
    static_cast<int>(prediction(prediction_row, q + h)) - 1;
  const int training_order =
    static_cast<int>(training(train_row, q + h)) - 1;
  double order_distance = 0.0;
  for (int coordinate = 0; coordinate < mapping.dimension; ++coordinate) {
    const double delta =
      mapping_value(mapping, prediction_order, coordinate) -
      mapping_value(mapping, training_order, coordinate);
    order_distance += delta * delta;
  }
  return std::exp(
    -parameters[q + h] * difference * difference - order_distance
  );
}

} // namespace

// [[Rcpp::export]]
NumericMatrix cpp_magp_mapping_matrix(
    const NumericVector& parameters,
    int q,
    int mapping_type) {
  MappingInfo mapping = build_mapping(parameters, q, mapping_type);
  NumericMatrix result(q, mapping.dimension);
  for (int level = 0; level < q; ++level) {
    for (int coordinate = 0; coordinate < mapping.dimension; ++coordinate) {
      result(level, coordinate) = mapping_value(mapping, level, coordinate);
    }
  }
  return result;
}

// [[Rcpp::export]]
NumericMatrix cpp_magp_covariance(
    const NumericMatrix& x,
    const NumericVector& parameters,
    int q,
    double tau,
    int mapping_type) {
  const int n = x.nrow();
  MappingInfo mapping = build_mapping(parameters, q, mapping_type);
  NumericMatrix covariance(n, n);
  double diagonal = tau;
  for (int h = 0; h < q; ++h) {
    diagonal += parameters[h];
  }
  for (int i = 0; i < n; ++i) {
    covariance(i, i) = diagonal;
    if ((i & 63) == 0) checkUserInterrupt();
    for (int j = i + 1; j < n; ++j) {
      double value = 0.0;
      for (int h = 0; h < q; ++h) {
        value += parameters[h] *
          component_base(x, i, j, h, q, parameters, mapping);
      }
      covariance(i, j) = value;
      covariance(j, i) = value;
    }
  }
  return covariance;
}

// [[Rcpp::export]]
NumericVector cpp_magp_gradient(
    const NumericMatrix& x,
    const NumericVector& parameters,
    int q,
    int mapping_type,
    const NumericMatrix& inverse_covariance,
    const NumericVector& alpha) {
  const int n = x.nrow();
  const int npar = parameters.size();
  MappingInfo mapping = build_mapping(parameters, q, mapping_type);
  NumericVector gradient(npar);

  // Only variance derivatives have nonzero covariance diagonals.
  double diagonal_weight = 0.0;
  for (int i = 0; i < n; ++i) {
    diagonal_weight += inverse_covariance(i, i) - alpha[i] * alpha[i];
  }
  for (int h = 0; h < q; ++h) {
    gradient[h] = diagonal_weight;
  }

  for (int i = 0; i < n; ++i) {
    if ((i & 63) == 0) checkUserInterrupt();
    for (int j = i + 1; j < n; ++j) {
      // The derivative matrix is symmetric, so combine both off-diagonal
      // contributions.
      const double weight =
        inverse_covariance(i, j) + inverse_covariance(j, i) -
        2.0 * alpha[i] * alpha[j];

      for (int h = 0; h < q; ++h) {
        const double difference = x(i, h) - x(j, h);
        const double base =
          component_base(x, i, j, h, q, parameters, mapping);
        const double covariance_component = parameters[h] * base;

        gradient[h] += weight * base;
        gradient[q + h] += weight *
          (-covariance_component * difference * difference);

        const int first_order = static_cast<int>(x(i, q + h)) - 1;
        const int second_order = static_cast<int>(x(j, q + h)) - 1;
        if (first_order == second_order) continue;

        for (int coordinate = 0;
             coordinate < mapping.dimension;
             ++coordinate) {
          const double first_value =
            mapping_value(mapping, first_order, coordinate);
          const double second_value =
            mapping_value(mapping, second_order, coordinate);

          const int first_parameter =
            mapping_parameter(mapping, first_order, coordinate);
          if (first_parameter >= 0) {
            gradient[2 * q + first_parameter] += weight *
              (-2.0 * covariance_component *
                (first_value - second_value));
          }

          const int second_parameter =
            mapping_parameter(mapping, second_order, coordinate);
          if (second_parameter >= 0) {
            gradient[2 * q + second_parameter] += weight *
              (-2.0 * covariance_component *
                (second_value - first_value));
          }
        }
      }
    }
  }
  return gradient;
}

// [[Rcpp::export]]
NumericMatrix cpp_magp_cross_covariance(
    const NumericMatrix& training,
    const NumericMatrix& prediction,
    const NumericVector& parameters,
    int q,
    int mapping_type) {
  const int n = training.nrow();
  const int n_new = prediction.nrow();
  MappingInfo mapping = build_mapping(parameters, q, mapping_type);
  NumericMatrix result(n_new, n);

  for (int i = 0; i < n_new; ++i) {
    if ((i & 255) == 0) checkUserInterrupt();
    for (int j = 0; j < n; ++j) {
      double value = 0.0;
      for (int h = 0; h < q; ++h) {
        value += parameters[h] * cross_component_base(
          training, prediction, j, i, h, q, parameters, mapping
        );
      }
      result(i, j) = value;
    }
  }
  return result;
}
