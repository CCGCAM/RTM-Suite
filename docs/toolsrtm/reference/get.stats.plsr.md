# Get Statistical Scores for PLSR Model

This function calculates statistical scores such as R-squared, RMSE,
MNMB,MB, FGE and MAE. for both training and testing datasets based on
predictions from a Partial Least Squares Regression (PLSR) model.

## Usage

``` r
get.stats.plsr(model, k, train, test, var)
```

## Arguments

- model:

  The trained ML model.

- k:

  The number of components used for prediction.

- train:

  The training dataset.

- test:

  The testing dataset.

- var:

  The variable of interest.

## Value

A data frame containing the statistical scores (R-squared, RMSE, MAE)
for both training and testing datasets.
