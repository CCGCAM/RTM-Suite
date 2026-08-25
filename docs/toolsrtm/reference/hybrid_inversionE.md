# Inversion of plant traits using ML models (Ensemble variant)

Inversion of plant traits using ML models (Ensemble variant)

## Usage

``` r
hybrid_inversionE(
  LUT = NULL,
  input = NULL,
  split = 0.8,
  setseed = NULL,
  collinearity = NULL,
  pattern = NULL,
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

- collinearity:

  collinearity-and-stepwise-vif-selection or CARS method implemmented,
  options='VIF' and 'CARS'

- pattern:

  Please indicate the number of bands with same pattern'B'.

- Field.data:

  dataframe with observations

- acron:

  acronynm for the observation measure: e.g., Cab_obsrv, where
  acron='\_observ' and Cab has same name as input

## Value

A list: `model` (the fitted ensemble model), `Stats` (accuracy
statistics on the held-out test split), and `Plot` (a ggplot scatter of
predicted vs. measured/simulated values for `input`).
