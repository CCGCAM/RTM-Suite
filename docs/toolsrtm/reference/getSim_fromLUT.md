# Simulate reflectance spectra while sweeping one trait, all others fixed

Builds a LUT where every trait is held at its default except `trait`,
which is swept from `nmin` to `nmax` in steps of `Interval`, runs the
canopy simulation, and plots the resulting spectra.

## Usage

``` r
getSim_fromLUT(
  trait = "Cab",
  nmin = 0,
  nmax = 100,
  Interval = 10,
  psoil = 0.5,
  model = "PROSAIL",
  method = "ggplot"
)
```

## Arguments

- trait:

  character. Name of the trait to sweep (e.g. `"Cab"`, `"LAI"`). Default
  `"Cab"`.

- nmin, nmax:

  numeric. Minimum/maximum value of `trait`'s sweep range.

- Interval:

  numeric. Step size between successive values of `trait` in the sweep.

- psoil:

  numeric (0-1). Soil brightness factor (0 = wet soil, 1 = dry soil),
  mixed from `ToolsRTM::dataSpec_PDB`'s dry/wet reference spectra.

- model:

  character. Canopy simulation to run: `"PROSAIL"`, `"INFORM"`, or
  `"PROSPECT"`.

- method:

  character. `"ggplot"` (default) or `"classical"` (base-R plotting) for
  the spectral plot.

## Value

If `method = "ggplot"`: a list with `LUT` (the swept-trait LUT used) and
`Plot` (the ggplot object). If `method = "classical"`: just the LUT
data.frame (the plot is drawn as a side effect via base graphics, not
returned).

## Examples

``` r
if (FALSE) { # \dontrun{
out <- getSim_fromLUT(trait = "Cab", nmin = 10, nmax = 80, Interval = 10, model = "PROSAIL")
out$Plot
} # }
```
