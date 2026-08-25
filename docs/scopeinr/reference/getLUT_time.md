# Get LUT table for SCOPE model

Get LUT table for SCOPE model

## Usage

``` r
getLUT_time(
  inputs = NULL,
  dataICOS = NULL,
  timeStart = "2018-01-01",
  timeEnd = "2018-01-30",
  freq.hour = "diurnal",
  nLUT = 100,
  SetFixedSeed = T
)
```

## Arguments

- inputs:

  a dataframe with the ranges and the distribution

- dataICOS:

  Flux network with radiation, air temperature ...

- timeStart:

  the initial time step in this format '2018-01-01'

- timeEnd:

  the initial time step in this format '2018-01-30'

- freq.hour:

  the diurnal time or hourly time options: 'diurnal', 'hourly' NULL at
  12.00 GTM

- nLUT:

  number of rows for each time step

- SetFixedSeed:

  Fixed seed

## Value

a dataframe with the LUT for Time series
