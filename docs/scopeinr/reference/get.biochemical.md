# get.biochemical

`get.biochemical` Calculates net assimilation rate A, fluorescence F
using biochemical model

This function calculates:

- stomatal resistance of a leaf or needle (s m-1)

- photosynthesis of a leaf or needle (umol m-2 s-1)

- fluorescence of a leaf or needle (fraction of fluor. in the dark)

## Usage

``` r
get.biochemical(data.leafbio, data.meteo, data.opts, fV, get.plots)
```

## Arguments

- data.leafbio:

  list. Leaf biochemical parameters: `Type` ('C3' or 'C4'),
  `stressfactor`, `Vcmax25`, `BallBerry0`, `BallBerrySlope`, `Rdparam`,
  `Kn0`/`Knalpha`/`Knbeta` (NPQ model parameters), optional `g_m`
  (mesophyll conductance, mol m^-2 s^-1 bar^-1), and (when temperature
  correction is enabled) the `TDP` sub-list of temperature-dependence
  parameters.

- data.meteo:

  list. Meteorological/leaf micro-environment characteristics: `Q`
  (absorbed PAR, umol photons m^-2 s^-1), `Cs` (CO2 concentration in the
  leaf boundary layer, ppm), `Temp` (leaf temperature, deg C or K), `eb`
  (vapour pressure in the leaf boundary layer, hPa), `Oa` (O2
  concentration, mmol/mol), `p` (air pressure, hPa).

- data.opts:

  data.frame. Simulation options; row 7 (`tempcor`) is a 0/1 flag for
  whether temperature correction of Vcmax and the rate constants is
  applied.

- fV:

  numeric. Scaling factor applied to `Vcmax25` (e.g. the canopy
  nitrogen/Vcmax profile factor estimated in the energy balance).

- get.plots:

  logical. If `TRUE`, returns diagnostic plots for stomatal conductance,
  assimilation, and Jmax/Vcmax rates.

## Value

A list `biochem_out` with (among others): `A` (net assimilation rate,
umol m^-2 s^-1), `Ci`/`Cc` (internal/chloroplast CO2 concentration,
ppm), `rcw` (stomatal resistance, s m^-1), `gs` (stomatal conductance),
`Vcmax`, `Rd`, `Ja` (actual electron transport rate), `ps`/`ps_rel`
(photochemical yield and degree of light saturation), `eta`
(fluorescence yield relative to dark-adapted), `qE`/`qQ`
(non-photochemical/photochemical quenching), `fs`/`fo`/`fm`
(steady-state/dark/light-saturated fluorescence yields),
`Kn`/`NPQ`/`Kf`/`Kp`/`Kd` (rate constants), and `SIF` (solar-induced
fluorescence, `fs * Q`).

## References

Farquhar et al. 1980, Collatz et al (1991, 1992), and: Dutta, D.,
Schimel, D. S., Sun, Y., Tol, C. V. D., & Frankenberg, C. (2019).
Optimal inverse estimation of ecosystem parameters from observations of
carbon and energy fluxes. Biogeosciences, 16(1), 77-103.

Van der Tol, C., Berry, J. A., Campbell, P. K. E., & Rascher, U. (2014).
Models of fluorescence and photosynthesis for interpreting measurements
of solar induced chlorophyll fluorescence. Journal of Geophysical
Research: Biogeosciences, 119(12), 2312-2327.

Bonan, G. B., Lawrence, P. J., Oleson, K. W., Levis, S., Jung, M.,
Reichstein, M., ... & Swenson, S. C. (2011). Improving canopy processes
in the Community Land Model version 4 (CLM4) using global flux fields
empirically inferred from FLUXNET data. Journal of Geophysical Research:
Biogeosciences, 116(G2).

## Author

    Joe Berry and Christiaan van der Tol, Ari Kornfeld (Original version in Matlab)

Carlos Camino (Ported version into R)

last updates Date: 21 Sep 2012.

Update: 20 Feb 2013. Update: Aug 2013: correction of L171: Ci \<-
Ci\*1e6/ p \* 1E3

Update: 2016-10 - (JAK) major rewrite to accomodate an iterative
solution to the Ball-Berry equation - also allows for g_m to be
specified for C3 plants, but only if Ci_input is provided.

Update: 25 Feb 2021: Temperature reponse functions by Dutta et al.
implemented

## Examples

``` r
if (FALSE) { # \dontrun{
biochem_out <- get.biochemical(data.leafbio, data.meteo, data.opts, fV = 1, get.plots = FALSE)
} # }
```
