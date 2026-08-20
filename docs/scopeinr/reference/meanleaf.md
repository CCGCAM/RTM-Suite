# Calculates the layer average and the canopy average of leaf properties per layer, per leaf angle and per leaf azimuth (36)

`meanleaf` averages a leaf-scale property, given as a 3D array over leaf
inclination, leaf azimuth and canopy layer, weighted by the leaf
inclination distribution function (`lidf`) and/or a per-layer weight
`Ps` (typically the fraction of sunlit leaf area).

## Usage

``` r
meanleaf(canopy, F_, canopy.choice, Ps)
```

## Arguments

- canopy:

  list. Canopy structural properties: `nlayers` (number of canopy
  layers), `nlincl` (number of leaf inclination classes), `nlazi`
  (number of leaf azimuth classes), `lidf` (leaf inclination
  distribution function weights).

- F\_:

  numeric array with dimensions (nli, nlazi, nl). Leaf-scale property to
  average, indexed by leaf inclination class, leaf azimuth class, and
  canopy layer.

- canopy.choice:

  character. Integration method: `'angles'` (integration over leaf
  angles only, one value per layer), `'layers'` (integration over layers
  only, weighted by `Ps`, single value), or `'angles_and_layers'`
  (integration over both leaf angles and layers, weighted by `lidf` and
  `Ps`).

- Ps:

  numeric vector of length nl. Fraction of sunlit leaf area per layer,
  used as the layer weight for
  `canopy.choice %in% c('layers', 'angles_and_layers')`.

## Value

numeric vector. `Fout_vertical`: for `canopy.choice = 'angles'`, one
averaged value per layer (length `nl`); for `'layers'`, a single
averaged value; for `'angles_and_layers'`, one averaged value per layer
(length `nl`).

## Details

Last update: 7 December 2007; update 11 February 2008 made modular
(Joris Timmermans); update 25 Feb 2013 Wout Verhoef: proposed name
change, removed globals and used canopy-structure for input.

## Author

    Christiaan van der Tol (Original version in Matlab)

Carlos Camino (Ported version into R)

## Examples

``` r
if (FALSE) { # \dontrun{
Fout <- meanleaf(canopy, F_, canopy.choice = 'angles_and_layers', Ps)
} # }
```
