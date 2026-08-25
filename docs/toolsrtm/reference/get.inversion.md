# get.inversion of the plant traits using different machine learning models

get.inversion of the plant traits using different machine learning
models

## Usage

``` r
get.inversion(
  data,
  depVar,
  inputs,
  algorithm = "PLSR",
  method.resampling = NULL,
  n.cores = NULL,
  seed = 123,
  n.samples = 500,
  save.model = FALSE,
  save.path = NULL,
  ...
)
```

## Arguments

- data:

  The data frame containing the variables.

- depVar:

  The dependent variable.

- inputs:

  The independent variables.

- algorithm:

  The type of machine learning model. Options: "PLSR" (Partial Least
  Squares Regression), "SVM" (Support Vector Machine), "RF" (Random
  Forest); "NN" (Neural Network); "GB" (Gradient boosting); "xGB"
  (eXtreme Gradient Boosting (XGBoost) with linear base learners);
  'Bayesian' ( Bayesian Additive Regression Trees); 'AdaBag' ( Bagged
  AdaBoost); "qLASSO" (Quantile Regression with LASSO penalty); "RVM"
  (Relevance Vector Machines (RVM) with linear kernel); 'BRNN' (Bayesian
  Regularized Neural Networks); "Ensemble" (Stacking Ensemble models)
  Default is "PLSR".

- method.resampling:

  The resampling method for controlling tht ML: Options are: "boot"
  (Bootstrapping); "boot632" (Bootstrapping-632); "optimism_boot";
  "boot_all"; "cv" (cross-Validation); "repeatedcv" (repeats k-fold
  cross-validation with 3 times); "LOOCV" (Leave-One-Out
  Cross-Validation with 3 times); "LGOCV"

- n.cores:

  The number of cores

- seed:

  The seed for reproducibility. Default is 123.

- n.samples:

  A integer with the number of samples used for tunning search
  (nsample/2) and create the ML model (n.sample)

- save.model:

  Logical indicating whether to save the trained models. Default is
  FALSE.

- save.path:

  Path to save the trained models. Required if save_models is TRUE.

- ...:

  Additional arguments (currently unused, reserved for future
  extensions).

## Value

A list containing predictions, statistics, and plots for the specified
machine learning model.

## Examples

``` r
if (FALSE) { # \dontrun{
get.inversion(data = my_data, depVar = "Cab", inputs = c("NDVI", "TCARI"), ML = "SVM", seed = 123)
} # }
```
