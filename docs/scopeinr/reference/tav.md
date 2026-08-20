# Stern's formula in Lekner & Dorf (1988) gives reflectance for alfa = 90 degrees

Stern's formula in Lekner & Dorf (1988) gives reflectance for alfa = 90
degrees

## Usage

``` r
tav(alfa, nr)
```

## Arguments

- alfa:

  numeric. Maximum incidence solid angle, in degrees.

- nr:

  numeric. Leaf refractive index.

## Value

angles: average transmissivity of a dielectric surface.

## Author

    Christiaan van der Tol (Original version in Matlab)

Carlos Camino (Ported version into R)

## Examples

``` r
tav(59, 1.4)
#> [1] 0.9641415
```
