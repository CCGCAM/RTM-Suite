# Spectral global sensitivity analysis (Sobol total index per wavelength)

Runs a canopy radiative transfer model many times while varying a set of
plant/soil traits, then computes the Sobol total sensitivity index
([`get.sobol.indices`](get.sobol.indices.md)) at each wavelength – i.e.
how much of the variance in simulated reflectance at that wavelength is
explained by each trait. Produces the data behind a stacked "Total SI
(percent) vs wavelength" plot (traits stacked to 100 percent at every
wavelength).

## Usage

``` r
get.spectral.sensitivity(
  n.samples = 1000,
  distribution = "Uniform",
  traits = c("N", "Cab", "EWT", "LMA", "LIDFa", "LAI"),
  leaf.model = "PROSPECT-D",
  canopy.model = "fourSAIL",
  rsoil.base = NULL,
  wl.step = 5,
  seed = 123,
  n.cores = NULL,
  chunk.size = 500,
  save.path = NULL
)
```

## Arguments

- n.samples:

  integer. Number of RTM simulations to run (also the total number of
  rows fed to [`get.sobol.indices`](get.sobol.indices.md), whose own `N`
  parameter is set to `n.samples / 2`, its internal split-half sample
  size). Default 1000, matching 500 used for the actual Sobol
  calculation.

- distribution:

  character. `"Uniform"` or `"Gaussian"` – which PDF each varied trait
  is drawn from (bounds/mean+sd come from `ToolsRTM::inputsPROSAIL`).

- traits:

  character vector. Which trait names (must exist in
  `ToolsRTM::inputsPROSAIL$variable`) to vary and attribute sensitivity
  to. Default: N, Cab, EWT, LMA, LIDFa, LAI – plus `SoilCoef` (a
  soil-brightness multiplier, 0.5-1.5, on the flat baseline soil
  spectrum) added automatically, matching the classic PROSAIL
  sensitivity figure (leaf structure, pigment, water, dry matter, leaf
  angle, LAI, soil).

- leaf.model, canopy.model:

  character. Passed to [`simulate_RTM`](simulate_RTM.md).

- rsoil.base:

  numeric. Baseline soil reflectance spectrum (before the `SoilCoef`
  multiplier); default flat 0.15 across 400-2500nm.

- wl.step:

  integer. Compute Sobol indices every `wl.step` nm instead of at every
  single nm, for speed (the Sobol calculation itself, not the RTM
  simulations, dominates runtime at full 1nm resolution). Default 5.

- seed:

  integer. Random seed for reproducibility.

- n.cores:

  integer. Number of cores for the parallel simulation loop (via
  `parallel`/`doParallel`, same pattern as `Scripts/Sensibility`).
  Default `parallel::detectCores() - 2`.

- chunk.size:

  integer. Simulations are run in chunks of this size (default 500)
  instead of all at once, matching the chunking pattern in
  `Scripts/1-getSCOPE-v3_withChunck.R`. This keeps memory bounded at
  large `n.samples` (5000, 20000, ...) and, if `save.path` is given,
  means a crash partway through only loses the in-progress chunk, not
  the whole run.

- save.path:

  character or `NULL`. If given, each chunk's raw simulated
  reflectance + LUT rows are saved to
  `file.path(save.path, "chunk_<i>.rds")` as soon as that chunk finishes
  (so results survive a crash/interrupt), and are read back from disk
  (not re-simulated) if the files already exist – re-running with the
  same `save.path` resumes instead of restarting.

## Value

A data.frame with columns `wavelength`, `trait`, and `STi_pct` (total
Sobol index, normalized to sum to 100% across traits at each wavelength)
– long format, ready for `ggplot2::geom_area(position = "stack")`.

## Examples

``` r
if (FALSE) { # \dontrun{
si_uniform <- get.spectral.sensitivity(n.samples = 1000, distribution = "Uniform")
si_normal  <- get.spectral.sensitivity(n.samples = 1000, distribution = "Gaussian")
} # }
```
