# Compute the numerical Jacobian of the Fluspect leaf model with respect to its parameters

Compute the numerical Jacobian of the Fluspect leaf model with respect
to its parameters

## Usage

``` r
get.numjacobian(x, spectral, inputLeaf, optipar)
```

## Arguments

- x:

  numeric vector. Leaf biochemistry/structure parameters at which to
  compute the numerical Jacobian.

- spectral:

  list. Spectral configuration used by the Fluspect model.

- inputLeaf:

  list. Leaf input parameters.

- optipar:

  list. Optical parameters used by the Fluspect model.

## Value

J: a numeric array of partial derivatives (reflectance/transmittance
with respect to each parameter in x).

## Author

    Wout Verhoef, Christiaan van der Tol, Joris Timmermans (Original version in Matlab)

Carlos Camino (Ported version into R)

## Examples

``` r
if (FALSE) { # \dontrun{
get.numjacobian(x, spectral, inputLeaf, optipar)
} # }
```
