# Inversion of plant traits using ML models

Inversion of plant traits using ML models

## Usage

``` r
hybrid_inversion(
  LUT = NULL,
  input = NULL,
  split = 0.8,
  setseed = NULL,
  method = NULL,
  collinearity = NULL,
  pattern = NULL,
  trans = NULL,
  Field.data = NULL,
  acron = NULL
)
```

## Arguments

- LUT:

  Dataset with inputs and Bands

- input:

  variable to estimate

- split:

  ratio between 0 and 1 for splitting the dataset in training and
  testing

- setseed:

  set random number

- method:

  Machine learning approach for estimating each plant traits, the
  options are 'SVM', 'RF' and 'LDA'

- collinearity:

  collinearity-and-stepwise-vif-selection or CARS method implemmented,
  options='VIF' and 'CARS'

- pattern:

  Please indicate the number of bands with same pattern'B'.

- trans:

  Please indicate is want a logarithm transformation to y variable .
  Default is T

- Field.data:

  dataframe with observations

- acron:

  acronynm for the observation measure: e.g., Cab_obsrv, where
  acron='\_observ' and Cab has same name as input

## Value

A list: `model` (the fitted ML model), `Stats` (accuracy statistics on
the held-out test split), and `Plot` (a ggplot scatter of predicted vs.
measured/simulated values for `input`). If `Field.data` is supplied, the
returned `Plot`/`Stats` are computed against the field observations
instead of the LUT's own test split.
