# this function produces Spatial maps of the spectral indices

this function produces Spatial maps of the spectral indices

## Usage

``` r
getSpectraIndices(
  rasterFiles = NULL,
  Sensor = "Sentinel2a",
  SpecEq = NULL,
  SpectraltoCompute = "All",
  factorR = NULL,
  path.export = NULL,
  single.bands = T
)
```

## Arguments

- rasterFiles:

  path with the images (in RasterBrick/RasterStack formats)

- Sensor:

  character with the sensor 'Sentinel2a'

- SpecEq:

  Character. Refer to an expression corresponding to the spectral index
  to compute

- SpectraltoCompute:

  List with the spectral index to compute. by default 'All'

- factorR:

  numeric. multiplying factor used to write reflectance in image
  (==10000 for S2)

- path.export:

  path_to save Bands

- single.bands:

  extract singles_bands, by deault is False

## Value

spectral indices
