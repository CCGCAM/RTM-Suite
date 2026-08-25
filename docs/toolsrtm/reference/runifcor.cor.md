# Generate a vector from a uniform distribution correlated with another distribution

Allow user to draw from a uniform distribution correlated with a user
specified vector

## Usage

``` r
runifcor.cor(x, rho)
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

Rough estimate

## Note

runif only works coming from uniform data

## References

modified from Eric Neuwirth, R-help mailing list post "Generating
uniformly distributed correlated data" (2011)

## Author

Jared E. Knowles

## Examples

``` r
x <- runif(1000)
y <- runifcor.cor(x, 0.2)
cor(x,y) # very close to 0.2
#> [1] 0.2152277
mean(y) 
#> [1] 0.5085712
sd(y)   
#> [1] 0.2936689
```
