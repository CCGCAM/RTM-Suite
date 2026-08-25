# get.high.temp.inhibtionC3 `get.high.temp.inhibtionC3` High Temperature Inhibition Function:The following function pertains to C3 photosynthesis

get.high.temp.inhibtionC3 `get.high.temp.inhibtionC3` High Temperature
Inhibition Function:The following function pertains to C3 photosynthesis

## Usage

``` r
get.high.temp.inhibtionC3(Tref, R, T, deltaS, deltaHd)
```

## Arguments

- Tref:

  numeric. Reference temperature, in Kelvin.

- R:

  numeric. Universal gas constant.

- T:

  numeric. Leaf temperature, in Kelvin.

- deltaS:

  numeric. Entropy term of the high-temperature inhibition function.

- deltaHd:

  numeric. Deactivation energy (enthalpy of high-temperature
  inhibition).

## Value

fHTv: high-temperature inhibition factor applied to C3 photosynthesis
rates.

## Examples

``` r
get.high.temp.inhibtionC3(Tref = 298.15, R = 8.314, T = 305, deltaS = 650, deltaHd = 200000)
#> [1] 0.720543
```
