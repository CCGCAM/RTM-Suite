# Run the combined B/Cx Fluspect leaf model from a packed parameter vector

Run the combined B/Cx Fluspect leaf model from a packed parameter vector

## Usage

``` r
get.calc_fluspect_bcar(params, spectral, leafbio, optipar)
```

## Arguments

- params:

  numeric vector of length 6, in order: Cab, Cdm, Cw, Cs, Cca, N - leaf
  biochemistry/structure parameters.

- spectral:

  list. Spectral configuration used by the Fluspect model.

- leafbio:

  list. Baseline leaf parameters; Cab/Cdm/Cw/Cs/Cca/N are overwritten
  from params before simulation.

- optipar:

  list. Optical parameters used by the Fluspect model.

## Value

A list with refl and tran (leaf reflectance/transmittance).

## Examples

``` r
if (FALSE) { # \dontrun{
get.calc_fluspect_bcar(params, spectral, leafbio, optipar)
} # }
```
