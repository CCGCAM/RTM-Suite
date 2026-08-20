# Build a canopy-model-ready soil reflectance spectrum from MARMIT

High-level wrapper around
[`get.marmit1`](get.marmit1.md)/[`get.marmit2`](get.marmit2.md) that
loads the measured dry-soil reference spectrum and water-optics
constants bundled with ToolsRTM, runs the MARMIT model, and resamples
the result onto a fixed wavelength grid so it can be dropped straight
into `rsoil` for [`foursail`](foursail.md), [`foursail2`](foursail2.md),
[`inform`](inform.md), or into [`SPART`](SPART.md)'s `rsoil` override.

## Usage

``` r
get.marmit.rsoil(
  database = "Bablet_2016",
  id = 1,
  version = "marmit1",
  L = 0.05,
  eps = 0.3,
  n_i = 1.53,
  k_i = 0.001,
  d_i = 5e-04,
  wl.out = 400:2500
)
```

## Arguments

- database:

  character. Soil database name (a folder under
  `extdata/marmit/databases/`). Default `"Bablet_2016"`, the only one
  bundled.

- id:

  integer. Soil ID within the database's index CSV (see the `ID`
  column). Default `1`.

- version:

  character. `"marmit1"` (default) or `"marmit2"`. MARMIT-2 additionally
  accounts for soil particle size/refractive index (`n_i`/`k_i`/`d_i`)
  and is generally more accurate for coarser soils; MARMIT-1 is simpler
  and matches the original 2018 paper.

- L:

  numeric. Thickness of the surface water layer, cm. Default `0.05`.

- eps:

  numeric (0-1). Fraction of the soil surface that is wet. Default
  `0.3`.

- n_i, k_i, d_i:

  numeric. MARMIT-2-only soil-particle parameters (real refractive
  index, imaginary refractive index, particle volume fraction). Ignored
  when `version = "marmit1"`. Defaults (`1.53`, `0.001`, `0.0005`) match
  the MARMIT Shiny app's defaults.

- wl.out:

  integer vector. Wavelength grid (nm) the output is resampled/padded
  onto. Default `400:2500` (matches [`foursail`](foursail.md)'s
  `spectrum.all = TRUE` / [`inform`](inform.md) grid); pass `400:2400`
  for Fluspect-leaf-model calls to
  [`foursail`](foursail.md)/[`foursail2`](foursail2.md).

## Value

A list:

- wavelength:

  the `wl.out` grid.

- rsoil.dry:

  dry-soil reflectance on `wl.out` (no MARMIT wetting applied – the raw
  measured reference).

- rsoil.wet:

  MARMIT-simulated reflectance on `wl.out`, ready to use as `rsoil`.

- SMC:

  estimated gravimetric soil moisture content (percent), from
  [`sigmoid.soil`](sigmoid.soil.md).

- params:

  a one-row data.frame recording `database`, `id`, `version`, `L`,
  `eps`, `n_i`, `k_i`, `d_i`.

## Details

Only the **Bablet 2016** soil database (Bablet et al., 2018) is bundled
with the package, to keep install size small. Other MARMIT databases
(Dupiau 2020, Humper 2015, Lesaignoux 2008, Liu 2002, Lobell 2002, Marcq
2012, Philpot 2014 – see
<https://pss-gitlab.math.univ-paris-diderot.fr/marmit/marmit>) can be
used by dropping a folder in the same layout (an index CSV
`<name>/<name>.csv` with columns `ID, Refl_file, SMCg, K, a, psi`, plus
`<name>/spectra/<Refl_file>` tab-separated `Wvl,R` files) under
`system.file("extdata", "marmit", "databases", package = "ToolsRTM")`
and passing `database = "<name>"`.

The dry-soil reference for a given `id` is the driest spectrum on file
for that soil (the row with the smallest `SMCg`), matching how the
original MARMIT Shiny app selects it. Wavelengths beyond the native
range of that spectrum (some databases stop at 2400 or 2490nm, not 2500)
are held constant at the last available value – the same boundary
behavior [`SPART`](SPART.md) uses for its own BSM soil beyond 2400nm.

## See also

[`get.marmit1`](get.marmit1.md), [`get.marmit2`](get.marmit2.md)

## Examples

``` r
if (FALSE) { # \dontrun{
soil <- get.marmit.rsoil(database = "Bablet_2016", id = 1, L = 0.05, eps = 0.3)
plot(soil$wavelength, soil$rsoil.wet, type = "l")
lines(soil$wavelength, soil$rsoil.dry, col = "grey50")

# Feed straight into fourSAIL (PROSPECT-D domain, full 400-2500nm)
LUT <- as.data.frame(getLUT(inputs = ToolsRTM::inputsPROSAIL, nLUT = 1, setseed = 1))
sim <- foursail(inputLUT = LUT, rsoil = soil$rsoil.wet, LeafModel = "PROSPECT-D")
} # }
```
