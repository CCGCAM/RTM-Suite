# Get first and total Sobol indices by Carlos Camino

This function calculates the first-order and total Sobol indices for
sensitivity analysis based on the input data and predictions.

## Usage

``` r
get.sobol.indices(data, output, N, normalize = FALSE)
```

## Arguments

- data:

  A data frame containing the input variables and the predicted output.

- output:

  A character string specifying the name of the column with the
  predictions.

- N:

  An integer specifying the number of samples to be used for Sobol
  sensitivity analysis.

- normalize:

  A boolean value; if TRUE, the Sobol indices will be normalized.

## Value

A data frame containing the first-order and total Sobol indices for each
input variable.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example usage
data <- data.frame(input1 = runif(100), input2 = runif(100), output = rnorm(100))
sobol_indices <- get.sobol.indices(data, output = "output", N = 1000, normalize = TRUE)
print(sobol_indices)
} # }
```
