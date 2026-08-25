# getPredicts

getPredicts

## Usage

``` r
getPredicts(
  model = NULL,
  type.model = "CNN",
  data = NULL,
  data.trans = NULL,
  data.Y = NULL,
  depVar = "Cab",
  scaler.depVar = NULL
)
```

## Arguments

- model:

  a keras model

- type.model:

  keras model type; the options are: 'CNN'; 'Hidden-layers'

- data:

  the dataset in matrix format for predicting depVar

- data.trans:

  PCA and normalized, options are avalaible

- data.Y:

  the Y variable: Avalaible options are: 'NULL' in real cases or a
  vector or matrix or data.frame

- depVar:

  variable name

- scaler.depVar:

  Scaler for Y

## Value

predictions
