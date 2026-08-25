# Generate a vector from a binomial distribution correlated with another distribution

Allow user to draw from a binomial distribution correlated with a user
specified vector

## Usage

``` r
rbinomcor(x, rho, scale = NULL, na.rm = FALSE)
```

## Arguments

- x:

  variable to draw from

- rho:

  correlation coefficient between x and result of function

- scale:

  a scale factor for the binomial draws

- na.rm:

  a logical indicating whether to fit the distribution excluding missing
  values or to fail on missing values

## Value

a vector of the same length as x drawn from a binomial distribution
correlated with x at the level of rho

## Details

Rough estimate

## Author

Jared E. Knowles
