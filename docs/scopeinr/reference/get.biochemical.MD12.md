# get.biochemical.MD12

`get.biochemical.MD12` Calculates:

- CO2 concentration in intercellular spaces (umol/mol == ppmv)

- leaf net photosynthesis (umol/m2/s) of C3 or C4 species

- fluorescence yield of a leaf (fraction of reference fluorescence yield
  in dark-adapted and un-stressed leaf)

Note: always use the prescribed units. Temperature can be either oC or K
Note: input can be single numbers, vectors, or n-dimensional matrices
Note: For consistency reasons, in C4 photosynthesis electron transport
rates under CO2-limited conditions are computed by inverting the
equation applied for light-limited conditions(Ubierna et al 2013). A
discontinuity would result when computing J from ATP requirements of Vp
and Vco, as a fixed electron transport partitioning is assumed for
light-limited conditions

Date: 21 Sep 2012 Update:

- 28 Jun 2013 Adaptation for use of Farquhar model of C3 photosynthesis
  (Farquhar et al 1980).

- 18 Jul 2013 Inclusion of von Caemmerer model of C4 photosynthesis (von
  Caemmerer 2000, 2013).

- 15 Aug 2013 Modified computation of CO2-limited electron transport in
  C4 species for consistency with light-limited value.

- 22 Oct 2013 Included effect of qLs on Jmax and electron transport
  value of kNPQs re-scaled in input as NPQs.

- 08 Jan 2019 (CvdT): minor modification to adjust to SCOPE_lite.

## Usage

``` r
get.biochemical.MD12(data.leafbio, data.meteo, fV, get.plots)
```

## Arguments

- data.leafbio:

  LUT table

- data.meteo:

  meteo characterisitics with L (Monin-Obukhov length), also carries the
  absorbed PAR (Q) and leaf temperature used internally.

- fV:

  fraction of Vcmax25 downregulated as function of cumulative absorbed
  PAR through the canopy (relative activity profile).

- get.plots:

  return plots for gs, assimilation and Jmax, Vcmax rate

## Value

the following parameters at eaf level:

- A in umol/m2/s which is the net assimilation rate of the leaves

- Ci in umol/mol which is the CO2 concentration in intercellular spaces
  (assumed to be the same as at carboxylation sites in C3 species)

- eta in (-) which is the amplification factor to be applied to PSII
  fluorescence yield spectrum relative to the dark-adapted, un-stressed
  yield calculated with either Fluspect or FluorMODleaf

## Author

Federico Magnani, with contributions from Christiaan van der Tol
(Original version in Matlab)

Carlos Camino (Ported version into R)

## Examples

``` r
if (FALSE) { # \dontrun{
out <- get.biochemical.MD12(data.leafbio, data.meteo, fV = 1, get.plots = FALSE)
} # }
```
