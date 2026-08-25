# Compute statistical metrics for model performance

This function evaluates training and testing performance using common
statistical accuracy metrics such as RMSE, MAE, MB, MNMB, FGE, and R².

## Usage

``` r
compute_stats(pred.train, train.obs, pred.test, test.obs)
```

## Arguments

- pred.train:

  Numeric vector of model predictions (training set)

- train.obs:

  Numeric vector of observed values (training set)

- pred.test:

  Numeric vector of model predictions (testing set)

- test.obs:

  Numeric vector of observed values (testing set)

## Value

A data.frame with Train/Test RMSE, MAE, MB, MNMB, FGE, and R².

## Examples

``` r
if (FALSE) { # \dontrun{
stats <- compute_stats(pred.train, train.obs, pred.test, test.obs)
} # }
```
