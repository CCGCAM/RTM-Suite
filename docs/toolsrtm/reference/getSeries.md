# get series in a dataframe from rasters

get series in a dataframe from rasters

## Usage

``` r
getSeries(
  pathRaster = NULL,
  shapefile = NULL,
  band_names = NULL,
  factorR = NULL,
  get.indices = T
)
```

## Arguments

- pathRaster:

  folder with the raster in tif format

- shapefile:

  shapefile file (.shp)

- band_names:

  names of the bands of the raster

- factorR:

  factor for the reflectance bands

- get.indices:

  A boolean is True, get spectral indices pre-define in getIndicesSE2,
  is not (FALSE) provide only bands

## Value

a dataframe with values
