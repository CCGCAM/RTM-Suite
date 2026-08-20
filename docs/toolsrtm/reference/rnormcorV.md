# Generate a vector from a normal distribution correlated with another distribution

Allow user to draw from a random normal distribution correlated with a
user specified vector

## Usage

``` r
rnormcorV(x, rho, ...)
```

## Arguments

- x:

  variable to draw from

- rho:

  correlation coefficient between x and result of function

## Value

a vector of the same length as x drawn from a normal distribution
correlated with x at the level of rho

## Details

Rough estimate, biased by known amount for now

## Author

Jared E. Knowles

## Examples

``` r
x <- rnorm(1000, 1, 1)
y <- rnormcorV(x, 0.2)
cor(x,y) # very close to 0.2
#> [1] 0.1349704
mean(y) # close to 0
#> [1] 0.1735345
sd(y)   # close to 1
#> [1] 1.01782
```
