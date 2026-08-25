# meanleaf.v2

`meanleaf.v2` averages a leaf-scale flux or property over leaf angle
classes and/or canopy layers, weighted by the leaf inclination
distribution function (`lidf`) and/or a layer weight `Ps` (typically the
sunlit/shaded fraction). This is a vectorized variant of
[`meanleaf`](meanleaf.md) used internally by
[`get.aggregator.ebal`](get.aggregator.ebal.md).

## Usage

``` r
meanleaf.v2(data.canopy, F_, canopy.choice, Ps)
```

## Arguments

- data.canopy:

  list. Canopy structural properties: `nlayers` (number of canopy
  layers), `nlincl` (number of leaf inclination classes), `nlazi`
  (number of leaf azimuth classes), `lidf` (leaf inclination
  distribution function weights).

- F\_:

  array or numeric vector. Leaf-scale flux/property to average; expected
  as an array with dimensions (nlincl, nlazi, nlayers) for
  `canopy.choice = 'angles'` or `'angles_and_layers'`, or as a numeric
  vector of length `nlayers` for `canopy.choice = 'layers'`.

- canopy.choice:

  character. Averaging mode: `'angles'` (average over leaf
  inclination/azimuth only), `'layers'` (weighted average over layers
  using `Ps`), or `'angles_and_layers'` (average over both).

- Ps:

  numeric vector of length nlayers, or scalar. Layer weight (e.g.
  fraction of sunlit leaf area) used when `canopy.choice` includes
  `'layers'`.

## Value

numeric. Single averaged value `Fout_vertical`.

## Author

    Christiaan van der Tol (Original version in Matlab)

Carlos Camino (Ported version into R)

## Examples

``` r
if (FALSE) { # \dontrun{
Fout <- meanleaf.v2(data.canopy, F_, canopy.choice = 'layers', Ps)
} # }
```
