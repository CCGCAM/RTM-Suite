# Satellite sensor band, SMAC atmospheric-correction and SRF metadata

`LANDSAT4.TM`, `LANDSAT5.TM`, `LANDSAT7.ETM`, `LANDSAT8.OLI`,
`Sentinel2A.MSI`, `Sentinel2B.MSI`, `Sentinel3A.OLCI`, `Sentinel3B.OLCI`
and `TerraAqua.MODIS` each bundle the per-band metadata, SMAC
atmospheric correction coefficients and spectral response function (SRF)
tables needed by
[`get.spectral.convolution()`](get.spectral.convolution.md) and the
sensor-convolution vignettes to simulate at-sensor reflectance for that
mission.

## Usage

``` r
LANDSAT4.TM

LANDSAT5.TM

LANDSAT7.ETM

LANDSAT8.OLI

Sentinel2A.MSI

Sentinel2B.MSI

Sentinel3A.OLCI

Sentinel3B.OLCI

TerraAqua.MODIS
```

## Format

A named list with entries that include `mission`/`name` (sensor
identifiers), `band_id_all` (band names), `res_spatials` (spatial
resolution per band), `rang_wvls` (nominal spectral range per band),
`swath_widths`, `revisit_days`/`revisit_time`, `band_width`,
`center_wvl` (band centre wavelengths), `SMAC_coef` (SMAC
atmospheric-correction coefficient table), `wl_smac`/`p_srf`/`wl_srf`
(or the `wvl_srf`/`p_srf_smac`/`wl_srf_smac` equivalents for the
Sentinel-2 entries) giving the per-band spectral response function
sampled on the SMAC/native wavelength grids, and `id_smac_in_all`
mapping SMAC bands to the full band list. Exact component names vary
slightly by sensor family; see the source list names for the specific
object.

An object of class `list` of length 17.

An object of class `list` of length 17.

An object of class `list` of length 17.

An object of class `list` of length 17.

An object of class `list` of length 17.

An object of class `list` of length 17.

An object of class `list` of length 17.

An object of class `list` of length 17.
