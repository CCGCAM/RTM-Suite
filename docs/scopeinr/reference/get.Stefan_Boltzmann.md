# Stefan-Boltzmann equation

`get.Stefan_Boltzmann` computes the total radiant exitance of a
blackbody at a given temperature, using the Stefan-Boltzmann law (\\H =
\sigma T^4\\).

## Usage

``` r
get.Stefan_Boltzmann(T_C)
```

## Arguments

- T_C:

  numeric vector. Temperature(s) of the emitting surface, in degrees
  Celsius.

## Value

numeric vector. Blackbody radiant exitance `H` (W m^-2) at each
temperature in `T_C`.

## Author

Carlos Camino

## Examples

``` r
if (FALSE) { # \dontrun{
H <- get.Stefan_Boltzmann(T_C = 20)
} # }
```
