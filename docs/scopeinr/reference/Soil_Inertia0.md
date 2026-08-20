# Calculate the soil thermal inertia from known soil thermal properties

`Soil_Inertia0` computes the soil thermal inertia `GAM` directly from
the soil specific heat, bulk density and thermal conductivity, used when
`options.soil_heat_method == 0` in the SCOPE energy balance.

## Usage

``` r
Soil_Inertia0(cs, rhos, lambdas)
```

## Arguments

- cs:

  numeric. Soil specific heat capacity (J kg^-1 K^-1).

- rhos:

  numeric. Soil bulk density (kg m^-3).

- lambdas:

  numeric. Soil thermal conductivity (W m^-1 K^-1).

## Value

numeric. Soil thermal inertia `GAM` (J m^-2 K^-1 s^-1/2), used to
compute the soil heat flux.

## Author

    Christiaan van der Tol (Original version in Matlab)

Carlos Camino (Ported version into R)

## Examples

``` r
if (FALSE) { # \dontrun{
GAM <- Soil_Inertia0(cs, rhos, lambdas)
} # }
```
