# Plot Predictions from PLSR Model

This function plots predictions from a Partial Least Squares Regression
(PLSR) model on both training and testing datasets.

## Usage

``` r
get.plot.plsr(model, k, data.train, data.test, var)
```

## Arguments

- model:

  The trained PLSR model.

- k:

  The number of components to use for prediction.

- data.train:

  The training dataset, must include the variable of interest.

- data.test:

  The testing dataset, must include the variable of interest.

- var:

  The variable of interest to plot (as a string).

## Value

A ggplot object containing the plot of predictions.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example usage
library(ggplot2)
# Assuming 'model', 'data.train', and 'data.test' are defined and 'var' is the target variable
plot <- get.plot.plsr(model, k = 2, data.train, data.test, var = "target_variable")
print(plot)

} # }
```
