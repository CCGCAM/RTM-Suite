# Aggregates MODTRAN data over SCOPE bands by averaging

`aggreg` reads a MODTRAN `.atm` output file and averages six relevant
atmospheric radiative transfer quantities (extraterrestrial solar
irradiance, hemispherical reflectance, direct and diffuse
transmittances, and path radiance) onto the (possibly coarser,
multi-resolution) SCOPE spectral band definition.

## Usage

``` r
aggreg(atmfile, SCOPEspec)
```

## Arguments

- atmfile:

  character. Path to a MODTRAN `.atm` text file with a header row,
  wavelength in column 2, and the relevant transmittance/radiance
  functions in columns 3-20.

- SCOPEspec:

  list. SCOPE spectral band definition, with elements `nreg` (number of
  spectral regions), `start` (start wavelength of each region, nm),
  `end` (end wavelength of each region, nm), and `res` (spectral
  resolution of each region, nm).

## Value

numeric matrix with `nwS` rows (total number of SCOPE bands across all
regions) and 6 columns. The six columns correspond to: (1)
`Eso*cos(tts)/pi`, (2) `rdd`, (3) `tss`, (4) `tsd`, (5) `tssrdd`, (6)
`La` (path radiance).

## Author

Carlos Camino

## Examples

``` r
if (FALSE) { # \dontrun{
M <- aggreg(atmfile, SCOPEspec)
} # }
```
