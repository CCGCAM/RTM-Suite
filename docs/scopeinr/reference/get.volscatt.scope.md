# get.volscatt.scope version 2.0 from SCOPE model

`volscatt.scope` calculates the scattering phase functions using a
radiative transfer model based on the optical properties of leaves.
Specifically, the model computes the fraction of radiation scattered in
the forward and backward directions for a given set of input parameters

## Usage

``` r
get.volscatt.scope(tts, tto, psi, ttli)
```

## Arguments

- tts:

  Sun: zenith angle in degrees

- tto:

  observation:zenith angle in degrees

- psi:

  Difference of azimuth angle between solar and viewing position

- ttli:

  leaf inclination array

## Value

The function returns a list of four items: "chi_s", "chi_o", "frho", and
"ftau". These items represent the scattering phase functions for direct
and diffuse radiation, respectively.

## Details

Date: 11 February 2008.

## Author

    Wout Verhoef, Joris Timmermans (Original version in Matlab)

Carlos Camino (Ported version into R)

## Examples

``` r
get.volscatt.scope(tts = 30, tto = 10, psi = 0, ttli = seq(5, 85, by = 10))
#> $chi_s
#> [1] 0.8627299 0.8365163 0.7848856 0.7094065 0.6123724 0.4967318 0.3891676
#> [8] 0.3412130 0.3207464
#> 
#> $chi_o
#> [1] 0.9810603 0.9512512 0.8925389 0.8067073 0.6963642 0.5648625 0.4161977
#> [8] 0.2548870 0.1239844
#> 
#> $frho
#> [1] 0.26951928 0.25421667 0.22545717 0.18670960 0.14264748 0.09858537 0.06081478
#> [8] 0.03304870 0.01583893
#> 
#> $ftau
#> [1] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [6] 0.000000e+00 9.769816e-04 1.970399e-03 6.324214e-05
#> 
```
