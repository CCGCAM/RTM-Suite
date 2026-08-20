# Soil thermal inertia method by Murray and Verhoef

`Soil_Inertia1` estimates the soil thermal inertia `GAM` from the soil
moisture content, using the empirical model of Murray and Verhoef (2007)
for thermal conductivity and a linear mixing model for volumetric heat
capacity. Used when `options.soil_heat_method == 1` in the SCOPE energy
balance.

## Usage

``` r
Soil_Inertia1(SMC)
```

## Arguments

- SMC:

  numeric. Soil moisture content (m^3 m^-3).

## Value

numeric. Soil thermal inertia `GAM` (J m^-2 K^-1 s^-1/2), used to
compute the soil heat flux.

## Author

    Christiaan van der Tol; Murray and Verhoef (Original version in Matlab)

Carlos Camino (Ported version into R)

## Examples

``` r
if (FALSE) { # \dontrun{
GAM <- Soil_Inertia1(SMC)
} # }
```
