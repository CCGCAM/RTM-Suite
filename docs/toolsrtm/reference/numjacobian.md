# num-jacobian function

num-jacobian function

## Usage

``` r
numjacobian(x, spectral, inputLeaf, optipar)
```

## Arguments

- x:

  numeric vector. Leaf biochemistry/structure parameters at which to
  compute the numerical Jacobian.

- spectral:

  list. Spectral configuration (wavelength ranges/resolution) used by
  the Fluspect model.

- inputLeaf:

  list. Leaf input parameters passed through to `calc_fluspect_bcar`.

- optipar:

  list. Optical parameters (refractive index, specific absorption
  coefficients) used by the Fluspect model.

## Value

numjacobian values

## Author

Wout Verhoef, Christiaan van der Tol, Joris Timmermans,

Ported to R: Carlos. Camino
