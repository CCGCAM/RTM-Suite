# LIBERTY model

Leaf radiative transfer model designed for conifer needles.

## Usage

``` r
liberty(inputLUT)
```

## Value

Length 4 list containing modeled reflectance, transmittance,
wavelengths, and R (?).

## Details

Note that these coefficients are based on the original values from
Dawson et al. (1998). Improved specific absorption coefficients are
available from later work by Di Vittorio (2009).

## References

Dawson, T. P., Curran, P. J., & Plummer, S. E. (1998). LIBERTY—Modeling
the Effects of Leaf Biochemical Concentration on Reflectance Spectra.
Remote Sensing of Environment, 65(1), 50–60.
https://doi.org/10.1016/S0034-4257(98)00007-8

Di Vittorio, A. V. (2009). Enhancing a leaf radiative transfer model to
estimate concentrations and in vivo specific absorption coefficients of
total carotenoids and chlorophylls a and b from single-needle
reflectance and transmittance. Remote Sensing of Environment, 113(9),
1948–1966. https://doi.org/10.1016/j.rse.2009.05.002

## Author

Alexey Shiklomanov and modifed by Carlos Camino
