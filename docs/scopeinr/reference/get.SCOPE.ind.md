# get simulations based on SCOPE model

Main module for calling the SCOPE model

## Usage

``` r
get.SCOPE.ind(
  LUT,
  options.SCOPE = data.opts,
  optipar,
  path.out,
  leaf.model = "fluspect-CX",
  canopy.model = "fourSAIL",
  get.outputs,
  get.plots = T
)
```

## Arguments

- LUT:

  leaf, biochemistry, viewing angles,meteo, and canopy properties needed
  for running SCOPE model.

- options.SCOPE:

  simulation options controlling the SCOPE run, such as the atmospheric
  correction and whether to correct the total irradiance if a specific
  value is provided instead of the usual Modtran output.

- optipar:

  leaf optical parameters used by the selected leaf model (e.g. Fluspect
  refractive index and specific absorption coefficients).

- path.out:

  folder for saving the SCOPE outputs.

- leaf.model:

  Not currently functional – see [`get.SCOPE`](get.SCOPE.md). Passing
  anything other than `'fluspect-CX'` triggers a warning.

- canopy.model:

  Not currently functional – see [`get.SCOPE`](get.SCOPE.md) (same
  underlying limitation: SCOPE's canopy engine isn't swappable yet).
  Passing anything other than `'fourSAIL'` triggers a warning.

- get.outputs:

  if get.outputs = 'ALL' all variables were retrieved; if get.outputs =
  'Main' the main variables were retrieved. By default SCOPE uses
  get.outputs = 'ALL', for processiong huge LUT is recommended
  get.outputs = 'Main'.

- get.plots:

  is true plot the intermediate plots

## References

Yang, P., E. Prikaziuk, W. Verhoef, and C. van der Tol. 2020. “SCOPE
2.0: A Model to Simulate Vegetated Land Surface Fluxes and Satellite
Signals.” Geoscientific Model Development Discussions 2020: 1–26.
https://doi.org/10.5194/gmd-2020-251.

Van der Tol, C., W. Verhoef, J Timmermans, A Verhoef, and Z Su. 2009.
“An Integrated Model of Soil-Canopy Spectral Radiances, Photosynthesis,
Fluorescence, Temperature and Energy Balance.” Biogeosciences 6 (12):
3109–29. https://doi.org/10.5194/bg-6-3109-2009.

## Author

    Wout Verhoef, Prikaziuk, Christiaan van der Tol, Joris Timmermans (Original version in Matlab)

Carlos Camino (Ported version into R)

## Examples

``` r
if (FALSE) { # \dontrun{
out <- get.SCOPE.ind(LUT, options.SCOPE = data.opts, optipar, path.out,
                      leaf.model = "fluspect-CX", canopy.model = "fourSAIL",
                      get.outputs = "Main", get.plots = FALSE)
} # }
```
