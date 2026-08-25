# extract spectral indices at Sentinel-2 resolution

extract spectral indices at Sentinel-2 resolution

## Usage

``` r
getIndicesSE2(df, sensor = "Sentinel-2a", df.data = NULL, fast.process = NULL)
```

## Arguments

- df:

  a dataframe with reflectance where each rows correspond with an
  spectrum

- sensor:

  Sensor options: 'Sentinel-2a', or 'Sentinel-2b'

- df.data:

  dataset with IDs that corresponde with each spectrum, is null is also
  enable

- fast.process:

  when the bands are ordered for SE2, please use fast.process = T,
  otherwise use False or nothing

## Value

a dataframe with indices and your dataset
