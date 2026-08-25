# Sentinel-2A/2B MSI measured spectral response functions

Full measured per-band spectral response function (SRF) curves for the
Sentinel-2A and Sentinel-2B MultiSpectral Instrument (MSI), sampled on a
common wavelength grid, as distributed by ESA.

## Usage

``` r
srf.sentinel2a

srf.sentinel2b
```

## Format

A data frame with 2301 rows (wavelength grid, `SR_WL`) and 14 columns:
`SR_WL` followed by one column per MSI band (B1-B12, B8A) giving that
band's measured relative spectral response at each wavelength.

An object of class `data.frame` with 2301 rows and 14 columns.
