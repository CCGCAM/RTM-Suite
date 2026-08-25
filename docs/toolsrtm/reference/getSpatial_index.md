# this function produces Spatial maps of the spectral index

this function produces Spatial maps of the spectral index

## Usage

``` r
getSpatial_index(
  rasterFiles = NULL,
  Sensor = "Sentinel2a",
  SpectraltoCompute = NULL,
  factorR = NULL
)
```

## Arguments

- rasterFiles:

  path with the image (in RasterBrick/RasterStack formats)

- Sensor:

  character with the sensor 'Sentinel2a'

- SpectraltoCompute:

  List with the spectral index to compute. by default 'All'

- factorR:

  numeric. multiplying factor used to write reflectance in image
  (==10000 for S2)

## Value

a spatial map.
