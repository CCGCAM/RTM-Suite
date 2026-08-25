# function for Spatial mapping of the main traits

function for Spatial mapping of the main traits

## Usage

``` r
getSpatialTrait(
  rasterFiles = NULL,
  ForestLayer = NULL,
  Sensor = NULL,
  saveFile = NULL,
  proj = 3035,
  shapeLayer = NULL,
  model.ML = NULL,
  trait = "Cab",
  factorR = NULL
)
```

## Arguments

- rasterFiles:

  path with the image (in RasterBrick/RasterStack formats) or a single
  stackRaster

- ForestLayer:

  path with the images of the forest mask. This code is adapted to from
  COPERNICUS Land Monitoring Service

- Sensor:

  by default is 'Sentinel2a', is not needed

- saveFile:

  path with the path to save the spatial trait mapping

- proj:

  projection of the files (shapefile and raster should be in same
  projection)

- shapeLayer:

  path with the shapefile with your study area. This shape is use for
  cropping images

- model.ML:

  Machine learning model

- trait:

  trais to estimate (actually for Cab and LAI)

- factorR:

  a numeric scaling factor applied to the raster reflectance values
  before index computation (e.g. to convert to 0-1 reflectance)

## Value

a list with: i) a scatter-plot between trait and best indicator, and ii)
spatial trait mappping mask with the forest map
