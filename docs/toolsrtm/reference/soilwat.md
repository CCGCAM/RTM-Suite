# soilwat function `soilwat` In this model it is assumed that the water film area is built up

Date: Version 1.0, September 2012 Wout Verhoef Update: Version 1.1, Jan
2020, Peiqi Yang

## Usage

``` r
soilwat(rdry, nw, kw, SMp, SMC, film)
```

## Arguments

- rdry:

  dry soil reflectance (NW,1)

- nw:

  refraction index of water in (NW,1)

- kw:

  absorption coefficient of water in (NW,1)

- SMp:

  soil moisture volume percentage (1,NS)

- SMC:

  soil moisture capacity (recommended 0.25) (1,1)

- film:

  effective optical thickness (deleff) of single water film (1,1)
  (recommended 0.015)

## Value

soil spectrum

## Author

    Wout Verhoef, Peiqi Yang, Christiaan van der Tol (Original version in Matlab)

Carlos Camino (Ported version into R)

## Examples

``` r
if (FALSE) { # \dontrun{
soilwat(rdry, nw, kw, SMp, SMC, film)
} # }
```
