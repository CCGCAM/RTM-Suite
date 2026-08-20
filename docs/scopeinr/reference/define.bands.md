# Define spectral regions used by SCOPEinR

`define.bands` defines the wavelength grid and spectral sub-regions used
throughout SCOPEinR (optical, thermal, PAR, excitation-fluorescence
matrix, etc.), following the SCOPE v1.40 spectral band definition (three
regions of increasing resolution: 400-2400 nm at 1 nm, 2500-15000 nm at
100 nm, and 16000-50000 nm at 1000 nm).

## Usage

``` r
define.bands()
```

## Value

A list `spectral` with (among others): `wlS` (full wavelength vector,
nm), `wlP` (PROSPECT range), `wlE` (excitation wavelengths for the E-F
fluorescence matrix), `wlF` (chlorophyll fluorescence emission
wavelengths), `wlO` (optical part), `wlT` (thermal part), `wlPAR` (PAR
range, 400-700 nm), and `SCOPEspec` (region boundaries and resolutions
used by [`aggreg`](aggreg.md) to read MODTRAN data).

## Author

    Wout Verhoef (Original version in Matlab)

Carlos Camino (Ported version into R)

## Examples

``` r
spectral <- define.bands()
```
