# Package index

## Fit and predict

Fit either mapping model and predict outcomes with uncertainty.

- [`magp2d_fit()`](https://twskr.github.io/magp/reference/magp2d_fit.md)
  : Fit a MaGP model with a two-dimensional sequence map
- [`magpfull_fit()`](https://twskr.github.io/magp/reference/magpfull_fit.md)
  : Fit a MaGP model with a full sequence map
- [`predict(`*`<magp2d>`*`)`](https://twskr.github.io/magp/reference/predict.magp.md)
  [`predict(`*`<magpfull>`*`)`](https://twskr.github.io/magp/reference/predict.magp.md)
  : Predict outcomes from a fitted MaGP model
- [`magp2d_rmse()`](https://twskr.github.io/magp/reference/magp2d_rmse.md)
  : Calculate root-mean-squared prediction error

## Construct an initial design

Build and assess quantitative-sequence starting designs.

- [`magp_initial_design()`](https://twskr.github.io/magp/reference/magp_initial_design.md)
  : Construct a quantitative-sequence initial design
- [`magp_sequence_design()`](https://twskr.github.io/magp/reference/magp_sequence_design.md)
  : Construct the sequence portion of an initial design
- [`magp_quantitative_design()`](https://twskr.github.io/magp/reference/magp_quantitative_design.md)
  : Construct the quantitative portion of an initial design
- [`magp_sequence_criterion()`](https://twskr.github.io/magp/reference/magp_sequence_criterion.md)
  : Evaluate a sequence initial design
- [`magp_quantitative_criterion()`](https://twskr.github.io/magp/reference/magp_quantitative_criterion.md)
  : Evaluate a quantitative Latin hypercube
- [`magp_joint_criterion()`](https://twskr.github.io/magp/reference/magp_joint_criterion.md)
  : Evaluate a complete quantitative-sequence initial design

## Choose later experiments

Use expected improvement for sequential optimization.

- [`magp_expected_improvement()`](https://twskr.github.io/magp/reference/magp_expected_improvement.md)
  : Calculate expected improvement for a fitted MaGP model
- [`magp_next_point()`](https://twskr.github.io/magp/reference/magp_next_point.md)
  : Find the next quantitative-sequence experiment
- [`magp_bayes_optimize()`](https://twskr.github.io/magp/reference/magp_bayes_optimize.md)
  : Run sequential Bayesian optimization with a MaGP surrogate
- [`magp_bayes_optimize_from_scratch()`](https://twskr.github.io/magp/reference/magp_bayes_optimize_from_scratch.md)
  : Start Bayesian optimization from a generated initial design

## Object summaries

Print concise summaries of fitted models, designs, and searches.

- [`print(`*`<magp2d>`*`)`](https://twskr.github.io/magp/reference/print.magp2d.md)
  : Summarize a fitted two-dimensional MaGP model
- [`print(`*`<magpfull>`*`)`](https://twskr.github.io/magp/reference/print.magpfull.md)
  : Summarize a fitted full-mapping MaGP model
- [`print(`*`<magp_initial_design>`*`)`](https://twskr.github.io/magp/reference/print.magp_initial_design.md)
  : Print a quantitative-sequence initial design
- [`print(`*`<magp_quantitative_design>`*`)`](https://twskr.github.io/magp/reference/print.magp_quantitative_design.md)
  : Print a quantitative initial design
- [`print(`*`<magp_sequence_design>`*`)`](https://twskr.github.io/magp/reference/print.magp_sequence_design.md)
  : Print a sequence initial design
- [`print(`*`<magp_next_point>`*`)`](https://twskr.github.io/magp/reference/print.magp_next_point.md)
  : Print a MaGP acquisition-search result
- [`print(`*`<magp_bayes_opt>`*`)`](https://twskr.github.io/magp/reference/print.magp_bayes_opt.md)
  : Print a MaGP Bayesian optimization result
