# get.resistances

`get.resistances` calculates the aerodynamic and boundary-layer
resistances between the soil, the canopy and the reference (measurement)
height, following the two-layer resistance scheme of Wallace and Verhoef
(2000), including a stability correction for non-neutral atmospheric
conditions based on the Monin-Obukhov length.

## Usage

``` r
get.resistances(data.soil, data.canopy, data.meteo)
```

## Arguments

- data.soil:

  list. Soil properties; must contain `rbs` (soil boundary-layer
  resistance, s m^-1).

- data.canopy:

  list. Canopy structural properties: `Cd` (leaf drag coefficient),
  `LAI` (leaf area index, m^2 m^-2), `rwc` (within-canopy aerodynamic
  resistance, s m^-1), `zo` (roughness length for momentum, m), `d`
  (zero-plane displacement height, m), `hc` (vegetation height, m),
  `leafwidth` (characteristic leaf width, m).

- data.meteo:

  list. Meteorological data: `z` (measurement height, m), `u` (wind
  speed at `z`, m s^-1), `L` (Monin-Obukhov length, m; see
  [`get.Monin.Obukhov`](get.Monin.Obukhov.md)).

## Value

A list `resist_out` with (among others): `ustar` (friction velocity, m
s^-1), `uz0` (wind speed at `z0m`, m s^-1), `Kh` (eddy diffusivity for
heat, m^2 s^-1), `rai`/`rar`/`rac` (aerodynamic resistance in the
inertial sublayer / roughness sublayer / canopy layer, s m^-1), `rws`
(aerodynamic resistance within the canopy down to the soil, s m^-1),
`raa` (total aerodynamic resistance above the canopy, s m^-1), `rawc`
(total resistance within the canopy air space, s m^-1), `raws` (total
resistance between the canopy air space and the soil, s m^-1).

## References

Wallace and Verhoef (2000) 'Modelling interactions in mixed-plant
communities: light, water and carbon dioxide', in: BruceMarshall, Jeremy
A. Roberts (ed), 'Leaf Development and Canopy Growth', Sheffield
Academic Press, UK. ISBN 0849397693

ustar: Tennekes, H. (1973) 'The logaritmic wind profile', J.Atmospheric
Science, 30, 234-238

Psih: Paulson, C.A. (1970), The mathematical representation of wind
speed and temperature in the unstable atmospheric surface layer. J.
Applied Meteorol. 9,857-861

## Author

Anne Verhoef, Christiaan van der Tol, Joris Timmermans (Original version
in Matlab)

Carlos Camino (Ported version into R)

Last update: 01 Feb 2008

## Examples

``` r
if (FALSE) { # \dontrun{
resist_out <- get.resistances(data.soil, data.canopy, data.meteo)
} # }
```
