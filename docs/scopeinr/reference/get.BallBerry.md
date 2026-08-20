# get.BallBerry `get.BallBerry` get BerryBall value from Berry Model

get.BallBerry `get.BallBerry` get BerryBall value from Berry Model

## Usage

``` r
get.BallBerry(Cs, RH, A, BallBerrySlope, BallBerry0, minCi, Ci_input)
```

## Arguments

- Cs:

  is CO2 at leaf surface

- RH:

  is relative humidity

- A:

  is Net assimilation in 'same units of CO2 as Cs'/m2/s

- BallBerrySlope:

  parameter for BallBerry model

- BallBerry0:

  parameter for BallBerry model

- minCi:

  minimum Ci as a fraction of Cs (in case RH is very low?)

- Ci_input:

  will use only for Ci == Ci_input

## Value

A list with gs (stomatal conductance) and Ci (intercellular CO2
concentration).

## Author

    Christiaan van der Tol (Original version in Matlab)

Carlos Camino (Ported version into R)

## Examples

``` r
get.BallBerry(Cs = 400, RH = 0.7, A = 15, BallBerrySlope = 9, BallBerry0 = 0.01, minCi = 0.3)
#> $gs
#> [1] 0.24625
#> 
#> $Ci
#> [1] 302.5381
#> 
```
