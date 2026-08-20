# Perform CARS-PLS analysis.

This function performs Competitive Adaptive Reweighted Sampling (CARS)
for partial least squares (PLS) analysis.

## Usage

``` r
get.cars.pls(
  X,
  y,
  nLV = 2,
  fold = 10,
  scale.pretreat = 1,
  iteration = 50,
  PartitionType = "interleaved"
)
```

## Arguments

- X:

  The sample matrix, where samples are in rows and variables are in
  columns.

- y:

  The response variable.

- nLV:

  The number of latent variables in PLS. Default is 2.

- fold:

  The number of segments for cross-validation. Default is 10.

- scale.pretreat:

  Whether to scale the variables. 1 for scaling, 0 for no scaling (only
  centered). Default is 1.

- iteration:

  The number of Monte Carlo samplings in CARS. Default is 50.

- PartitionType:

  The partition type for cross-validation: "random", "consecutive", or
  "interleaved". Default is "interleaved".

## Value

A list containing the results of the CARSPLS analysis.

## Author

Ported by Carlos Camino; orignal code oin matlab by Yizeng Liang, and
Hongdong Li

## Examples

``` r
if (FALSE) { # \dontrun{
# X: sample matrix (rows = samples, columns = variables); y: response vector
get.cars.pls(X, y)
get.cars.pls(X, y, nLV = 3, fold = 5, scale.pretreat = 0)
} # }
```
