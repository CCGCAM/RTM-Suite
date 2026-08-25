# Get Statistical Scores between measurd and predicted single variable

This function calculates statistical scores such as R-squared, RMSE,
MNMB,MB, FGE and MAE.

## Usage

``` r
get.stats.ind(df, depVar = NULL, depVar.pred = NULL)
```

## Arguments

- df:

  a dataframe; the dataset.

- depVar:

  A character; the measured variable name

- depVar.pred:

  A character; the predicted variable name

## Value

A data frame containing the statistical scores (R-squared, RMSE, MAE)
for both training and testing datasets.
