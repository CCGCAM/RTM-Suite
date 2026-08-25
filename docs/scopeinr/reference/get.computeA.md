# get.computeA `get.computeA` Compute the net CO2 assimilation rate using the Farquhar model Note: even though computeA() is written as a separate function, the code is, in fact, executed exactly this point in the file (i.e. between the previous if clause and the next section

get.computeA `get.computeA` Compute the net CO2 assimilation rate using
the Farquhar model Note: even though computeA() is written as a separate
function, the code is, in fact, executed exactly this point in the file
(i.e. between the previous if clause and the next section

## Usage

``` r
get.computeA(
  Ci,
  Type,
  g_m,
  Vs_C3,
  MM_consts,
  Rd,
  Vcmax,
  Gamma_star,
  Je,
  effcon,
  atheta,
  kpepcase
)
```

## Arguments

- Ci:

  is internal CO2 concentration

- Type:

  is photosynthetic pathway type ("C3" or "C4")

- g_m:

  is mesophyll conductance

- Vs_C3:

  is maximum carboxylation rate for C3 plants

- MM_consts:

  is Michaelis-Menten constants for Rubisco and oxygenase reactions

- Rd:

  is dark respiration rate

- Vcmax:

  is maximum carboxylation rate

- Gamma_star:

  CO2 compensation point in the absence of mitochondrial respiration

- Je:

  electron transport rate in the chloroplasts

- effcon:

  electron transport to carboxylation efficiency

- atheta:

  curvature parameter for the response of electron transport to
  irradiance

- kpepcase:

  the relative limitation of electron transport by irradiance

## Value

A: net CO2 assimilation rate (and related intermediate variables,
depending on Type/g_m branch taken).

## Details

Even though `computeA()`-style logic is written as a separate function
here, the calculation is executed exactly at this point in the original
SCOPE model's control flow (i.e., between the previous if clause and the
next section) — this R port keeps that same ordering.
