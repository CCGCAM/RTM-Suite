# get.temperature.functionC3 `get.temperature.functionC3` Temperature Correction Functions:The following function pertains to C3 photosynthesis

get.temperature.functionC3 `get.temperature.functionC3` Temperature
Correction Functions:The following function pertains to C3
photosynthesis

## Usage

``` r
get.temperature.functionC3(Tref, R, Temp, deltaHa)
```

## Arguments

- Tref:

  numeric. Reference temperature, in Kelvin.

- R:

  numeric. Universal gas constant.

- Temp:

  numeric. Leaf temperature, in Kelvin.

- deltaHa:

  numeric. Activation energy of the temperature-dependent process.

## Value

fTv: Arrhenius-type temperature correction factor.

## Examples

``` r
get.temperature.functionC3(Tref = 298.15, R = 8.314, Temp = 305, deltaHa = 65330)
#> [1] 1.807444
```
