# `get.SCOPE.plots` get simulations based on SCOPE model

`get.SCOPE.plots` get simulations based on SCOPE model

## Usage

``` r
get.SCOPE.plots(path.files, plant.trait, get.plots)
```

## Arguments

- path.files:

  folder with main tables from the SCOPE simulation.

- plant.trait:

  plant trait

- get.plots:

  type plot for comparing, options are: 'fluorescence', 'reflectance',
  'radiance'

## Author

Carlos Camino (Ported version into R)

## Examples

``` r
if (FALSE) { # \dontrun{
get.SCOPE.plots(path.files = "outs", plant.trait = "Cab", get.plots = "reflectance")
} # }
```
