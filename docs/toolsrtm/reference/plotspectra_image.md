# Plot the spectra in points over the image

Plot the spectra in points over the image

## Usage

``` r
plotspectra_image(
  Raster = NULL,
  n_spectra = 1,
  sensor = "Sentinel2a",
  factor = NULL,
  method = "ggplot"
)
```

## Arguments

- Raster:

  a brick or a Stack

- n_spectra:

  number of spectra for selecting in the image

- sensor:

  the sensor, for the moment, this function works for Sentinel-2 images

- factor:

  the factor for reflectance

- method:

  or ggplot or classical method, by default is taken 'ggplot'

## Value

a dataframe with the selected point and also show a plot with the
seoctral signal of each poin
