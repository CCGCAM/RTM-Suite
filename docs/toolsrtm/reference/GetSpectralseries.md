# Extract all bands from SE images (adapted to NETCDF dataset)

Extract all bands from SE images (adapted to NETCDF dataset)

## Usage

``` r
GetSpectralseries(
  netCDFs = NULL,
  bands = NULL,
  shapefile = NULL,
  factorSE = 1/10000,
  Indices = T
)
```

## Arguments

- netCDFs:

  List of NetCDF

- bands:

  a vector with the names of the inputs of the NetCDF

- shapefile:

  shapefile (point)

- factorSE:

  factor to apply in reflectance files (point)

- Indices:

  if TRUE or null, the function estimate the Indices for each date.

## Value

spectral indices
