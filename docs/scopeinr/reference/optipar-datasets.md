# Leaf optical parameters (PROSPECT/Fluspect, various parameterizations)

Wavelength-dependent leaf-level optical coefficients (refractive index,
specific absorption coefficients per pigment, PSI/PSII quantum yield
spectra) required by the Fluspect/PROSPECT leaf models. Several versions
are bundled, corresponding to different published parameterizations.

## Usage

``` r
optipar

optipar.2015

optipar2017.ProspectD

optipar2020.prospectD.BSM2019

optipar2021.Pro.CX
```

## Format

A data frame/list with one row per wavelength and one column per optical
coefficient.

An object of class `data.frame` with 2001 rows and 15 columns.

An object of class `list` of length 15.

An object of class `list` of length 16.

An object of class `list` of length 18.

## Source

`optipar`, `optipar.2015`: Feret et al. 2008 (PROSPECT-4);
`optipar2017.ProspectD`: Feret et al. 2017 (PROSPECT-D);
`optipar2020.prospectD.BSM2019`: Feret et al. 2017 + BSM soil model
(Yang et al. 2019); `optipar2021.Pro.CX`: PROSPECT-PRO (Feret et al.
2021).
