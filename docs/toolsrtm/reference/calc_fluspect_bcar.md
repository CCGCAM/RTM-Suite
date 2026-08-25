# calc_fluspect_bcar function

calc_fluspect_bcar function

## Usage

``` r
calc_fluspect_bcar(params, spectral, leafbio, optipar)
```

## Arguments

- params:

  numeric vector of length 6, in order: Cab (chlorophyll content), Cdm
  (dry matter content), Cw (water content), Cs (senescent material/brown
  pigments), Cca (carotenoid content), N (leaf structure parameter).

- spectral:

  list. Spectral configuration (wavelength ranges/resolution) used by
  the Fluspect model.

- leafbio:

  list. Leaf biochemistry/structure parameters; `Cab`, `Cdm`, `Cw`,
  `Cs`, `Cca`, and `N` are overwritten from `params` before simulation.

- optipar:

  list. Optical parameters (refractive index, specific absorption
  coefficients) used by the Fluspect model.

## Value

rfl and trans
