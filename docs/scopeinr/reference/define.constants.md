# Define physical constants used by SCOPEinR

`define.constants` returns a list of the physical constants used
throughout the SCOPEinR model (radiative transfer, energy balance and
biochemical sub-models), analogous to the
[`SCOPEinR::constants`](constants.md) table.

## Usage

``` r
define.constants()
```

## Value

A list with named numeric elements: `A` (Avogadro's number, mol^-1), `h`
(Planck's constant, J s), `c` (speed of light, m s^-1), `cp` (specific
heat of dry air, J kg^-1 K^-1), `R` (molar gas constant, J mol^-1 K^-1),
`rhoa` (specific mass of air, kg m^-3), `g` (gravitational acceleration,
m s^-2), `kappa` (Von Karman constant, dimensionless), `MH2O` (molecular
mass of water, g mol^-1), `Mair` (molecular mass of dry air, g mol^-1),
`MCO2` (molecular mass of CO2, g mol^-1), `sigmaSB` (Stefan-Boltzmann
constant, W m^-2 K^-4), `deg2rad` (degrees-to-radians conversion
factor), `C2K` (Celsius-to-Kelvin offset, K).

## Author

    Wout Verhoef (Original version in Matlab)

Carlos Camino (Ported version into R)

## Examples

``` r
const <- define.constants()
```
