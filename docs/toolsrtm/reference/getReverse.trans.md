# Reverse a caret preProcess transformation and round the result

Reverse a caret preProcess transformation and round the result

## Usage

``` r
getReverse.trans(preProc = NULL, data = NULL, digits = 10)
```

## Arguments

- preProc:

  a "preProcess" object created by
  [`caret::preProcess()`](https://rdrr.io/pkg/caret/man/preProcess.html),
  describing the transformation to reverse.

- data:

  a matrix. Transformed data to convert back to its original scale.

- digits:

  integer. Number of decimal places to round the reversed values to.
  Default 10.

## Value

reverse trans
