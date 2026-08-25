# Angers leaf biochemistry/optics reference database

A bundled subset of the Angers leaf-optical-properties experimental
database (leaf biochemical traits paired with measured reflectance and
transmittance spectra), used as example/validation data for leaf model
inversion.

## Usage

``` r
LeafDB.Angers
```

## Format

A named list with components `DataBioch` (data frame of measured leaf
biochemical traits, one row per sample), `lambda` (wavelength vector for
the spectra), `Refl`/`Tran` (measured reflectance and transmittance
spectra, one row per sample) and `nbSamples` (sample count).
