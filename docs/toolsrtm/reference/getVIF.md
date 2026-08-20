# Calculate Variance Inflation Factor (VIF)

This function calculates the Variance Inflation Factor (VIF) for the
predictor variables in a linear regression model to assess
multicollinearity.

## Usage

``` r
getVIF(in_frame, thresh = 10, trace = T, ...)
```

## Arguments

- in_frame:

  A data frame containing the predictor variables. The dependent
  variable should not be included in this frame.

- thresh:

  A numeric value indicating the threshold for VIF; predictors with VIF
  greater than this threshold will be flagged as having
  multicollinearity. Default is 10.

- trace:

  A boolean value; if TRUE, the function will print information about
  the VIF calculations and any predictors that exceed the threshold.

- ...:

  Additional arguments passed to other methods (not used in this
  function).

## Value

A data frame containing the VIF values for each predictor variable.

## Examples

``` r
# Example data frame
df <- data.frame(x1 = rnorm(100), x2 = rnorm(100), x3 = rnorm(100))
df$x2 <- df$x1 + rnorm(100, sd = 0.1)  # Introduce multicollinearity
vif_results <- getVIF(in_frame = df, thresh = 5, trace = TRUE)
#> Registered S3 methods overwritten by 'fmsb':
#>   method    from
#>   print.roc pROC
#>   plot.roc  pROC
#>  var vif             
#>  x1  93.1930644924128
#>  x2  93.5379967655569
#>  x3  1.03613553381576
#> 
#> removed:  x2 93.538 
#> 
print(vif_results)
#> [1] "x1" "x3"
# The following VIF function were extracted from 
#https://beckmw.wordpress.com/2013/02/05/collinearity-and-stepwise-vif-selection/
```
