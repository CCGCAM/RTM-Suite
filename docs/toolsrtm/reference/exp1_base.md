# Compute Exponential Integral E1(x)

This function computes the exponential integral E1(x) for each element
in the input vector.

## Usage

``` r
exp1_base(x)
```

## Arguments

- x:

  A numeric vector of values.

## Value

A numeric vector containing the result of the exponential integral E1(x)
for each value in `x`.

## Examples

``` r
exp1_base(1:10)  # Compute E1 for values 1 to 10
#>  [1] 2.193839e-01 4.890051e-02 1.304838e-02 3.779331e-03 1.148288e-03
#>  [6] 3.600799e-04 1.154808e-04 3.766531e-05 1.244725e-05 4.156932e-06
```
