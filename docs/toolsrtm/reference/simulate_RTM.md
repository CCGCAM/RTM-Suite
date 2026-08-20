# Simulate canopy reflectance with any supported leaf model + any supported canopy model

A single entry point connecting ToolsRTM's leaf models to its canopy
models. There are exactly THREE independent canopy models in this
package: `"fourSAIL"`, `"foursail2"`, and `"INFORM"`. Everything else
you might see referenced internally (`foursail.inform`, `foursail.inf`,
`foursail_t_o`, `foursail_t_s`) are implementation details INFORM uses
internally — not standalone models you call directly. An earlier version
of this dispatcher incorrectly exposed those as top-level options;
that's fixed here.

## Usage

``` r
simulate_RTM(
  inputLUT,
  rsoil,
  leaf.model = "PROSPECT-PRO",
  canopy.model = "fourSAIL",
  PROSPECTversion = "PRO",
  LUT_GB = NULL,
  optipar = NULL,
  ...
)
```

## Arguments

- inputLUT:

  one row of a LUT (e.g. from [`getLUT()`](getLUT.md)).

- rsoil:

  soil reflectance spectrum, as required by the chosen canopy model.

- leaf.model:

  character. One of "PROSPECT-PRO", "PROSPECT-D", "Liberty",
  "Fluspect-B", "Fluspect-B-Cx" — all 5 work with all 3 canopy models.

- canopy.model:

  character. One of `"fourSAIL"`, `"foursail2"`, `"INFORM"` — the only
  three real canopy models in this package.

- PROSPECTversion:

  character. Legacy `foursail2`-only switch, superseded by `leaf.model`;
  kept for backward compatibility.

- LUT_GB:

  data.frame. Only used when `canopy.model = "foursail2"`: reference
  leaf parameters for the green and brown vegetation fractions (row 1 =
  green, row 2 = brown). If NULL, foursail2 uses its own built-in
  defaults.

- optipar:

  optical parameters, passed through when `leaf.model` is a Fluspect
  variant.

- ...:

  additional arguments passed through to the selected canopy model
  function.

## Value

A list with element `rsot` (bi-directional reflectance factor) always
present, plus `rdot`/`rddt`/`rsdt` when the selected canopy model
provides them (`"fourSAIL"`/`"foursail2"`: all four; `"INFORM"`: only
`rsot`, since INFORM's forest-level BRF has no separate components – the
others are `NULL`).

## Examples

``` r
if (FALSE) { # \dontrun{
inputs <- ToolsRTM::inputsPROSAIL
LUT <- as.data.frame(ToolsRTM::getLUT(inputs = inputs, nLUT = 1, setseed = 1234))
sim <- simulate_RTM(inputLUT = LUT[1, ], rsoil = rsoil, leaf.model = "PROSPECT-PRO", canopy.model = "fourSAIL")
} # }
```
