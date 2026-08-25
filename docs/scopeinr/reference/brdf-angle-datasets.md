# BRDF viewing/illumination angle configurations

Predefined sets of solar/viewing zenith and relative azimuth angles used
to compute the bidirectional reflectance distribution function (BRDF) of
the canopy. `brdf_angles_no_oversampling` omits the finer angular
oversampling used by `brdf_angles`/`brdf_angles2`.

## Usage

``` r
brdf_angles

brdf_angles2

brdf_angles_no_oversampling
```

## Format

A data frame with columns for solar zenith, viewing zenith and relative
azimuth angles (degrees).

An object of class `data.frame` with 333 rows and 2 columns.

An object of class `data.frame` with 55 rows and 2 columns.
