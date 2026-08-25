# EnMAP hyperspectral sensor band characteristics

Per-channel centre wavelength and full-width-half-maximum (FWHM) table
for the EnMAP hyperspectral imager, used by the sensor-convolution
functions to build a Gaussian spectral response for EnMAP bands.

## Usage

``` r
EnMap.characteristics
```

## Format

A data frame with 242 rows (one per EnMAP channel) and 4 columns:
`Sensor`, `channel` (band index), `center` (centre wavelength, nm) and
`fwhm` (full width at half maximum, nm).
