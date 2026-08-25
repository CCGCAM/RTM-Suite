# get.brdf

`get.brdf` simulates the bidirectional reflectance (and, optionally, the
directional thermal, fluorescence and xanthophyll-related radiance) of
the canopy over a large set of viewing angles, by repeatedly calling
[`getRTMo`](getRTMo.md) (and, depending on `data.opts`,
`get.RTMt.planck`, `get.RTMf`, `get.RTMz`) for each requested angle
combination. The viewing angles requested in `data.directional` are
extended with a set of angles for hot-spot oversampling and
principal-plane oversampling, and duplicate angle combinations are
removed before the loop.

## Usage

``` r
get.brdf(
  data.spectral,
  data.angles,
  data.rad,
  data.directional,
  atmo,
  data.soil,
  data.leafopt,
  data.leafbio,
  data.canopy,
  data.gap,
  data.meteo,
  data.thermal,
  data.bcu,
  data.bch,
  data.opts,
  get.plots = F
)
```

## Arguments

- data.spectral:

  list. Spectral band definitions used throughout the model (e.g. `wlS`,
  `wlF`, `wlT`).

- data.angles:

  list. Base observation geometry, with element `tts` (solar zenith
  angle, degrees) used to construct the hot-spot oversampling angles.

- data.rad:

  list. Radiation fluxes and top-of-canopy irradiance computed
  previously by `getRTMo`, passed through to the
  thermal/fluorescence/xanthophyll sub-models.

- data.directional:

  list. User-requested directional observation geometry, with elements
  `psi` (relative azimuth angles, degrees) and `tto` (viewing zenith
  angles, degrees).

- atmo:

  list or data.frame. Atmospheric input (MODTRAN transmittance functions
  or precomputed `Esun_`/`Esky_`), passed to `getRTMo`.

- data.soil:

  list. Soil reflectance and related soil properties.

- data.leafopt:

  list. Leaf optical properties (reflectance/transmittance spectra)
  computed by the leaf optical model.

- data.leafbio:

  list. Leaf biochemical/biophysical parameters (`Cab`, `Cw`, `Cdm`,
  `N`, ...).

- data.canopy:

  list. Canopy structural properties (LAI, leaf inclination
  distribution, hot-spot parameter, etc.).

- data.gap:

  list. Gap fraction probabilities (`Ps`, `Po`, `Pso`) needed to combine
  sunlit/shaded contributions in the thermal, fluorescence and
  xanthophyll sub-models.

- data.meteo:

  list. Meteorological forcing data (e.g. air temperature, vapour
  pressure).

- data.thermal:

  list. Component temperatures from the energy balance: `Tcu` (sunlit
  canopy), `Tch` (shaded canopy), `Tsu` (sunlit soil), `Tsh` (shaded
  soil), all in deg C.

- data.bcu:

  list. Biochemical outputs for sunlit leaves, including `eta` (relative
  fluorescence emission efficiency) and `Kn` (non-photochemical
  quenching rate constant).

- data.bch:

  list. Biochemical outputs for shaded leaves, with the same structure
  as `data.bcu`.

- data.opts:

  data.frame. Model configuration options; rows used here are
  `calc_planck` (compute thermal radiance spectrum), `calc_fluor`
  (compute chlorophyll fluorescence), and `calc_xanthophyllabs` (include
  xanthophyll de-epoxidation effect on reflectance).

- get.plots:

  logical. If `TRUE`, diagnostic plots are produced by the underlying
  RTM calls. Default `FALSE`.

## Value

A list `directional` with matrices indexed by wavelength (rows) and
viewing angle (columns), including:

- psi, tto:

  numeric vectors. The (deduplicated) relative azimuth and viewing
  zenith angles actually simulated.

- refl\_:

  TOC reflectance spectrum per angle.

- rso\_:

  Directional-directional BRDF per angle.

- Lo\_:

  TOC outgoing radiance per angle.

- LoF\_:

  TOC fluorescence radiance per angle (if `calc_fluor` is enabled).

- Lot\_:

  TOC thermal radiance per angle (if `calc_planck` is enabled).

- Eoutte, BrightnessT:

  Placeholder outputs (currently left at their initialized zero values).

## Author

Carlos Camino

## Examples

``` r
if (FALSE) { # \dontrun{
directional <- get.brdf(data.spectral, data.angles, data.rad, data.directional, atmo,
                         data.soil, data.leafopt, data.leafbio, data.canopy, data.gap,
                         data.meteo, data.thermal, data.bcu, data.bch, data.opts)
} # }
```
