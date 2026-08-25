# Fluspect specific absorption coefficients (2003/2008 parameterisation)

Leaf optical coefficient table used by the Fluspect
fluorescence-emission model implementation, including chlorophyll
fluorescence excitation efficiency spectra.

## Usage

``` r
optipar
```

## Format

A data frame with 2001 rows (400-2400 nm) and 20 columns including
wavelength (`wl`), refractive index (`nr`), specific absorption
coefficients (`Kab`, `Kca`, `Ks`, `Kw`, `Kdm`), fluorescence excitation
efficiencies (`phiI`, `phiII`), carotenoid violaxanthin/ zeaxanthin
coefficients (`KcaV`, `KcaZ`), anthocyanin coefficient (`Kant`) and
generalised soil/senescence vectors (`GSV.1`, `GSV.2`).
