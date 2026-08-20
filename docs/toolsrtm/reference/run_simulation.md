# Run a vegetation reflectance/energy-balance simulation, choosing models interactively

The single entry point tying ToolsRTM and SCOPEinR together. Two
simulation types are supported:

## Usage

``` r
run_simulation(
  simulation.type = "optical",
  inputLUT,
  rsoil = NULL,
  leaf.model = "PROSPECT-PRO",
  canopy.model = "fourSAIL",
  optipar = NULL,
  options.SCOPE = NULL,
  parallel = FALSE,
  get.outputs = "ALL",
  get.plots = FALSE,
  ...
)
```

## Arguments

- simulation.type:

  character. "optical" (ToolsRTM only, canopy reflectance) or "scope"
  (SCOPEinR, full energy balance). Default "optical".

- inputLUT:

  one row (for "optical") or full multi-row LUT (for "scope") of
  biophysical parameters.

- rsoil:

  soil reflectance spectrum. Required for "optical"; ignored for "scope"
  (SCOPEinR builds its own from `options.SCOPE`).

- leaf.model:

  character. Leaf model name. For "optical": "PROSPECT-PRO",
  "PROSPECT-D", "Liberty", "Fluspect-B", "Fluspect-B-Cx". For "scope":
  SCOPEinR's own naming, e.g. "fluspect-CX".

- canopy.model:

  character. Canopy model name. For "optical": "fourSAIL", "foursail2",
  "INFORM" (the only three real canopy models in ToolsRTM). For "scope":
  SCOPEinR's own naming, e.g. "fourSAIL".

- optipar:

  optical parameters. For "optical", passed to Fluspect variants if
  used. For "scope", passed to SCOPEinR (e.g.
  [`SCOPEinR::optipar2021.Pro.CX`](https://rdrr.io/pkg/SCOPEinR/man/optipar-datasets.html)).

- options.SCOPE:

  SCOPE run-options table (only used when `simulation.type = "scope"` –
  see SCOPEinR's own documentation for its structure).

- parallel:

  logical. For "scope" only: run in parallel via
  [`SCOPEinR::get.SCOPE.parallel()`](https://rdrr.io/pkg/SCOPEinR/man/get.SCOPE.parallel.html)
  instead of
  [`SCOPEinR::get.SCOPE()`](https://rdrr.io/pkg/SCOPEinR/man/get.SCOPE.html).
  Default FALSE.

- get.outputs:

  character. For "scope" only: which SCOPE outputs to return, e.g.
  "ALL". Default "ALL".

- get.plots:

  logical. For "scope" only: whether SCOPEinR should generate its own
  diagnostic plots. Default FALSE.

- ...:

  additional arguments passed through to the underlying simulate_RTM() /
  SCOPEinR call.

## Value

For "optical": the output of [`simulate_RTM()`](simulate_RTM.md)
(reflectance components). For "scope": the output of SCOPEinR's
`get.SCOPE()`/`get.SCOPE.parallel()` (reflectance, fluorescence,
radiance, and energy-balance fluxes).

## Details

**"optical"** – canopy reflectance only, no atmosphere/energy balance.
Runs entirely within ToolsRTM via [`simulate_RTM()`](simulate_RTM.md).
Fast, good for building large LUTs for inversion.

**"scope"** – full Soil-Canopy-Observation-Photochemistry-Energy balance
simulation (reflectance + fluorescence + radiance + energy fluxes), via
SCOPEinR. Requires the SCOPEinR package to be installed – this function
checks for it explicitly and gives an informative error (with the real
GitLab install command) rather than a cryptic "could not find function"
if it's missing, since SCOPEinR is a companion package, not a hard
dependency of ToolsRTM.

This function does not invent any simulation logic of its own: the
"optical" path calls the existing, tested
[`simulate_RTM()`](simulate_RTM.md); the "scope" path calls SCOPEinR's
own `get.SCOPE()` / `get.SCOPE.parallel()` using the same argument names
confirmed from real teaching scripts using both packages together.

## Examples

``` r
if (FALSE) { # \dontrun{
# Optical-only simulation (ToolsRTM alone, fast, no energy balance)
inputs <- ToolsRTM::inputsPROSAIL
LUT <- as.data.frame(ToolsRTM::getLUT(inputs = inputs, nLUT = 1, setseed = 1234))
sim <- run_simulation(
  simulation.type = "optical",
  inputLUT = LUT[1, ], rsoil = rsoil,
  leaf.model = "PROSPECT-PRO", canopy.model = "fourSAIL"
)

# Full SCOPE simulation with energy balance (needs SCOPEinR installed)
scope_opts <- read.table("Tables/inputs/setoptions.csv", header = TRUE, sep = ",")
sim_scope <- run_simulation(
  simulation.type = "scope",
  inputLUT = LUT, options.SCOPE = scope_opts,
  optipar = SCOPEinR::optipar2021.Pro.CX,
  leaf.model = "fluspect-CX", canopy.model = "fourSAIL",
  parallel = TRUE, get.outputs = "ALL"
)
} # }
```
