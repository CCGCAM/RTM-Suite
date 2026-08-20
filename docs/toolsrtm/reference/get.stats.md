# Get Statistical Scores for ML Models

This function calculates statistical scores such as R-squared, RMSE,
MNMB,MB, FGE and MAE. for both training and testing datasets based on
predictions from a ML model.

## Usage

``` r
get.stats(model, train, test, var)
```

## Arguments

- model:

  The trained ML model.

- train:

  The training dataset.

- test:

  The testing dataset.

- var:

  The variable of interest.

## Value

A data frame containing the statistical scores (R-squared, RMSE, MAE)
for both training and testing datasets.
