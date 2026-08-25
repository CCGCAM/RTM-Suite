# Multi-sensor spectral band centre/bounds lookup table

A combined lookup table of nominal band bounds and average centre
wavelength across the sensors supported by the package's sensor
convolution and index functions.

## Usage

``` r
sensor.characteristics
```

## Format

A data frame with 343 rows (one per sensor/band combination) and 5
columns: `Sensor`, `channel` (band name), `lb`/`ub` (lower/upper
wavelength bound of the band, nm) and `average` (band centre wavelength,
nm).
