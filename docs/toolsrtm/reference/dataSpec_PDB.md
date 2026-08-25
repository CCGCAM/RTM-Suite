# PROSPECT-family leaf optical specific absorption coefficients (PDB set)

Wavelength-indexed specific absorption/refractive-index coefficients
used by the PROSPECT-PRO / PROSPECT-D leaf optical model implementations
in this package (see [`prospect_DB()`](prospect_DB.md)).

## Usage

``` r
dataSpec_PDB
```

## Format

A data frame with 2101 rows (400-2500 nm, 1 nm steps) and 12 columns:
`wavelength`, `Refrac_leafm` (leaf refractive index), specific
absorption coefficients for chlorophyll (`SC_chl`), carotenoids
(`SC_car`), anthocyanins (`SC_Anth`), brown pigments (`SC_Brwon`), water
(`SC_Cw`) and dry matter (`SC_Cm`), plus direct/diffuse solar irradiance
and dry/wet soil reflectance reference spectra.
