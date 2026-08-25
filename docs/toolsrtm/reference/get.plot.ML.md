# This function plots predictions from a machine learning (ML) model on both training and testing datasets.

This function plots predictions from a machine learning (ML) model on
both training and testing datasets.

## Usage

``` r
get.plot.ML(model, data.train, data.test, var)
```

## Arguments

- model:

  The trained ML model.

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
plot <- get.plot.ML(model, data.train, data.test, var = "target_variable")
print(plot)

} # }
```
