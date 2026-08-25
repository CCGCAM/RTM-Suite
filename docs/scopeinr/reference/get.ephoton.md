# get.ephoton

`get.ephoton` calculates the energy content (J) of 1 photon of
wavelength lambda (m)

## Usage

``` r
get.ephoton(lambda, constants)
```

## Arguments

- lambda:

  numeric vector. Wavelength(s) of the photon(s) (m).

- constants:

  data.frame. Physical constants table (as used throughout SCOPEinR),
  must contain rows for `h` (Planck's constant, J s) and `c` (speed of
  light, m s^-1) in columns `constant`/`value`.

## Value

numeric vector. Energy content `E` of one photon at each wavelength
`lambda` (J).

## Author

     Wout Verhoef, Christiaan van der Tol, Joris Timmermans (Original version in Matlab)

Carlos Camino (Ported version into R)

## Examples

``` r
if (FALSE) { # \dontrun{
E <- get.ephoton(lambda, constants)
} # }
```
