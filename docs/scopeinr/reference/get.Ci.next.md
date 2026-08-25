# get.Ci.next `get.Ci.next` Function to calculate the difference between "guessed" Ci (Ci_in) and Ci computed using BB after computing A Test-function for iteration (note that it assigns A in the function's context.) As with the next section, this code can be read as if the function body executed at this point. (if iteration was used). In other words, A is assigned at this point in the file (when iterating).

get.Ci.next `get.Ci.next` Function to calculate the difference between
"guessed" Ci (Ci_in) and Ci computed using BB after computing A
Test-function for iteration (note that it assigns A in the function's
context.) As with the next section, this code can be read as if the
function body executed at this point. (if iteration was used). In other
words, A is assigned at this point in the file (when iterating).

## Usage

``` r
get.Ci.next(Ci_in, Cs, RH, minCi, BallBerrySlope, BallBerry0, A_fun, ppm2bar)
```

## Arguments

- Ci_in:

  numeric. The "guessed" intercellular CO2 concentration to test.

- Cs:

  numeric. CO2 concentration at the leaf surface.

- RH:

  numeric. Relative humidity (0-1).

- minCi:

  numeric. Minimum Ci as a fraction of Cs.

- BallBerrySlope:

  numeric. Slope parameter of the Ball-Berry stomatal conductance model.

- BallBerry0:

  numeric. Intercept (minimum conductance) parameter of the Ball-Berry
  model.

- A_fun:

  function. Assimilation function, called as `A_fun(Cs)` to compute net
  assimilation.

- ppm2bar:

  numeric. Conversion factor from ppm to bar for CO2 partial pressure.

## Value

A list with err (difference between guessed and Ball-Berry-computed Ci)
and Ci_out (the Ball-Berry result).

## Examples

``` r
if (FALSE) { # \dontrun{
get.Ci.next(Ci_in = 280, Cs = 400, RH = 0.7, minCi = 0.3, BallBerrySlope = 9,
            BallBerry0 = 0.01, A_fun = function(Cs) list(A = 15), ppm2bar = 1e-6)
} # }
```
