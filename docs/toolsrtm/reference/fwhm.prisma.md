# PRISMA hyperspectral sensor band FWHM table

Per-channel centre wavelength and full-width-half-maximum (FWHM) table
for the PRISMA hyperspectral imager, used to build a Gaussian spectral
response for PRISMA bands.

## Usage

``` r
fwhm.prisma
```

## Format

A data frame with 234 rows (one per PRISMA channel) and 3 columns: `QB`
(band quality/type flag), `wavelength` (centre wavelength, nm) and
`fwhm` (full width at half maximum, nm).
