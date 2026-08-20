# get.heatfluxes `get.heatfluxes` Calculates latent and sensible heat flux

get.heatfluxes `get.heatfluxes` Calculates latent and sensible heat flux

## Usage

``` r
get.heatfluxes(ra, rs, Tc, ea, Ta, e_to_q, Ca, Ci)
```

## Arguments

- ra:

  aerodynamic resistance for heat s m-1

- rs:

  stomatal resistance s m-1

- Tc:

  leaf temperature oC

- ea:

  vapour pressure above canopy hPa

- Ta:

  air temperature above canopy oC

- e_to_q:

  conv. from vapour pressure to abs hum hPa-1

- Ca:

  ambient CO2 concentration umol m-3

- Ci:

  intercellular CO2 concentration umol m-3

## Value

a list with lEc latent heat flux of a leaf W m-2

## Author

    Christiaan van der Tol (Original version in Matlab)

Carlos Camino (Ported version into R)

## Examples

``` r
if (FALSE) { # \dontrun{
out <- get.heatfluxes(ra, rs, Tc, ea, Ta, e_to_q, Ca, Ci)
} # }
```
