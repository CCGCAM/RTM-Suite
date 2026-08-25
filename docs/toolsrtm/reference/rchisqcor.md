# Generate a vector from a chi-square distribution correlated with another distribution

Allow user to draw from a random chi-square distribution correlated with
a user specified vector

## Usage

``` r
rchisqcor(x, rho, na.rm = FALSE)
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
x <- rnorm(1000, 1, 1)
y <- rchisqcor(x, 0.2)
cor(x,y) # very close to 0.2
#> [1] 0.2207573
mean(y) 
#> [1] 3.120634
sd(y)   
#> [1] 2.937076
```
