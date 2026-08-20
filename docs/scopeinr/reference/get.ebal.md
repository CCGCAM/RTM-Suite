# get.ebal `get.ebal` Calculates the energy balance of a vegetated surface

date 26 Nov 2007 (CvdT) updates:

- 29 Jan 2008 (JT & CvdT) converted into a function

- 11 Feb 2008 (JT & CvdT) improved soil heat flux and temperature
  calculation

- 14 Feb 2008 (JT) changed h in to hc (as h=Avogadro\`s constant)

- 31 Jul 2008 (CvdT) Included Pntot in output

- 19 Sep 2008 (CvdT) Converted F0 and F1 from units per aPAR into units
  per iPAR

- 07 Nov 2008 (CvdT) Changed layout

- 18 Sep 2012 (CvdT) Changed Oc, Cc, ec

- Feb 2012 (WV) introduced structures for variables

- Sep 2013 (JV, CvT) introduced additional biochemical model

- 10 Dec 2019 (CvdT) made a light version (layer averaged fluxes)

## Usage

``` r
get.ebal(
  data.rad,
  data.gap,
  data.meteo,
  data.soil,
  data.canopy,
  data.leafbio,
  data.leafopt,
  data.spectral,
  data.opts,
  integrate.layer,
  k.maxit,
  get.plots
)
```

## Arguments

- data.rad:

  incident radiation

- data.gap:

  probabilities of direct light penetration and viewing

- data.meteo:

  meteo characterisitics

- data.soil:

  soil properties

- data.canopy:

  canopy properties

- data.leafbio:

  leaf biochemical parameters

- data.leafopt:

  leaf biochemical parameters

- data.spectral:

  spectral information for the model

- data.opts:

  options for running the model

- integrate.layer:

  how estimate the sum of layer (angles, layers or angles_layers)

- k.maxit:

  maximum number for iteration on the balance model

- get.plots:

  is true plot the intermediate plots

## Value

a list with:

- iter numerical parameters used in the iteration for energy balance
  closure.

- fluxes energy balance, turbulent, and CO2 fluxes.

- rad radiation spectra.

- thermal temperatures, aerodynamic resistances and friction velocity.

- bcu, bch leaf biochemical outputs for sunlit and shaded leaves,

## Author

    Christiaan van der Tol,  Joris Timmermans (Original version in Matlab)

Carlos Camino (Ported version into R)

## Examples

``` r
if (FALSE) { # \dontrun{
out <- get.ebal(data.rad, data.gap, data.meteo, data.soil, data.canopy,
                 data.leafbio, data.leafopt, data.spectral, data.opts,
                 integrate.layer = "angles_layers", k.maxit = 100,
                 get.plots = FALSE)
} # }
```
