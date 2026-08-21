# Generate a vector from a Weibull distribution correlated with another distribution

Allow user to draw from a Weibull distribution correlated with a user
specified vector

## Usage

``` r
rweibullcor(x, rho, na.rm = FALSE)
```

## Arguments

- x:

  variable to draw from

- rho:

  correlation coefficient between x and result of function

- na.rm:

  a logical indicating whether to fit the distribution excluding missing
  values or to fail on missing values

## Value

a vector of the same length as x drawn from a normal distribution
correlated with x at the level of rho

## Details

Rough estimate

## Author

Jared E. Knowles

## Examples

``` r
x <- rnorm(1000)
y <- rweibullcor(x, 0.2)
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
cor(x,y) # very close to 0.2
#> [1] 0.1817967
mean(y) 
#> [1] 0.5097917
sd(y)  
#> [1] 0.334904
```
