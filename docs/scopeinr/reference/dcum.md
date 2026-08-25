# dcum function

`dcum` computes the cumulative leaf inclination distribution value at a
given angle `theta`, following the ellipsoidal LIDF model of Verhoef;
for `a >= 1` an analytical solution is used, otherwise the value is
found by fixed-point iteration.

## Usage

``` r
dcum(a, b, theta)
```

## Arguments

- a:

  numeric. LIDFa parameter, controls the average leaf slope.

- b:

  numeric. LIDFb parameter, controls the distribution's bimodality.

- theta:

  numeric. Leaf inclination angle (degrees).

## Value

numeric. Cumulative LIDF value `f` at angle `theta`, dimensionless, in
\[0,1\].

## Details

edit 2017-12-28: change sampling of angles to match with dcum.m

## Author

    Joris Timmermans, Christiaan van der Tol  (Original version in Matlab)

Carlos Camino (Ported version into R)

## Examples

``` r
if (FALSE) { # \dontrun{
f <- dcum(a = -0.35, b = -0.15, theta = 40)
} # }
```
