# get.aggregator.ebal

`get.aggregator.ebal` aggregates a leaf-scale flux, computed separately
for sunlit and shaded leaves, into a single LAI-scaled canopy-integrated
total flux, by internally calling [`meanleaf.v2`](meanleaf.v2.md) for
each contribution and combining them weighted by the fraction of sunlit
leaf area `Fs`.

## Usage

``` r
get.aggregator.ebal(
  LAI,
  sunlit_flux,
  shaded_flux,
  Fs,
  data.canopy,
  canopy.choice
)
```

## Arguments

- LAI:

  numeric. Total (one-sided) leaf area index of the canopy (m2 m-2).

- sunlit_flux:

  array with dimensions (nlincl, nlazi, nlayers) (or as required by
  `canopy.choice`). Leaf-scale flux for sunlit leaves.

- shaded_flux:

  numeric vector of length nlayers. Leaf-scale flux for shaded leaves
  (averaged over leaf angle, one value per layer).

- Fs:

  numeric vector of length nlayers. Fraction of sunlit leaf area at each
  canopy layer (`Ps`); the shaded contribution is weighted by `1 - Fs`.

- data.canopy:

  list. Canopy structural properties: `nlayers` (number of canopy
  layers), `nlincl` (number of leaf inclination classes), `nlazi`
  (number of leaf azimuth classes), `lidf` (leaf inclination
  distribution function weights).

- canopy.choice:

  character. Aggregation mode passed to [`meanleaf.v2`](meanleaf.v2.md):
  one of `'angles'`, `'layers'`, or `'angles_and_layers'`.

## Value

numeric. LAI-scaled canopy-integrated total flux `flux_tot`, combining
the sunlit and shaded contributions.

## Author

    Christiaan van der Tol (Original version in Matlab)

Carlos Camino (Ported version into R)

## Examples

``` r
if (FALSE) { # \dontrun{
flux_tot <- get.aggregator.ebal(LAI, sunlit_flux, shaded_flux, Fs,
                                 data.canopy, canopy.choice = 'angles_and_layers')
} # }
```
