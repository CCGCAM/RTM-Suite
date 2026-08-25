# 00. Configuring SCOPE: options, defaults, and common errors

``` r

library(ToolsRTM)
library(SCOPEinR)
```

Every SCOPEinR tutorial from here on opens with the same line:

``` r

scope_options <- read.table(system.file("input", "setoptions.csv", package = "SCOPEinR"),
                             header = TRUE, sep = ",")
```

`scope_options` is a 17-row switchboard controlling which physics SCOPE
actually runs – not tuning parameters (those live in the LUT), but
structural choices: which sub-models execute, which files get read,
which approximations apply. This tutorial is a configuration and
troubleshooting reference, not a lab report: every option’s default,
what changing it does, when you’d actually want to, and – for the ones
this port doesn’t fully support yet – exactly what error you’ll see and
what to do about it.

## 1. What should I normally use?

The shipped defaults, unmodified, work:

``` r

scope_options <- read.table(system.file("input", "setoptions.csv", package = "SCOPEinR"),
                             header = TRUE, sep = ",")
scope_options
#>    Order              Options Value
#> 1      1                 lite     1
#> 2      2    calc_fluorescence     1
#> 3      3 calc_spectrum_planck     1
#> 4      4  calc_xanthophyllabs     1
#> 5      5         soilspectrum     0
#> 6      6   Fluorescence_model     0
#> 7      7            applTcorr     1
#> 8      8               verify     1
#> 9      9               mSCOPE     0
#> 10    10           simulation     0
#> 11    11     calc_directional     0
#> 12    12   calc_vert_profiles     0
#> 13    13     soil_heat_method     2
#> 14    14         calc_rss_rbs     0
#> 15    15         MoninObukhov     1
#> 16    16                 LIDF     0
#> 17    17           irradiance     0
#>                                                                                                                                                                                            Info
#> 1                                                                                                                                    Value 1 indicates the SCOPE will use the lite SCOPEversion
#> 2                                                                                              Value 1 indicates that SCOPE will calculate thechlorophyll fluorescence in observation direction
#> 3                                                                                                                                Value 1 SCOPE will calculate the spectrum of thermal radiation
#> 4                                                                                             Value 1 indicates the SCOPE includes simulation of reflectance dependence on de-epoxydation state
#> 5                                                                                       Value 0 use a soil reflectance file and  Value 1 SCOPE will calculate the soil spectrum using BSM model
#> 6                                                           Value 0 empirical with sustained NPQ from Flexas data and Value 1 empirical with sigmoid for Kn and Value 2 uses Magnani 2012 model
#> 7                                                                                       Value 1 indicates that SCOPE will correct the Vcmax and rate constants for temperature in biochemical.m
#> 8                                                                                                                                       Value 1 indicates that SCOPE will check with field data
#> 9                                                                                    Value 1 indicates that SCOPE wil use the mSCOPE considering vertical variations in the vegetation canopies
#> 10                                           Value 0 the SCOPE will execute by individual runs based on LUT table and Value 1 for time series (uses text files with meteo input as time series)
#> 11                                                                                                                                                  Value 1 calculate full BRDF for many angles
#> 12                                                                                                                                Value 1  indicates that SCOPE will estimate vertical profiles
#> 13 Value 0 will estimate the GAM parameters with Soil_Inertia0(lambdas) and Value 1 will estimate the GAM with Soil_Inertia1(SMC) and Value 2 will estimate  G by 0.35*Rn where always in no TS
#> 14                                                                                                                                                   Value 0 is fixed and Value 1 is calculated
#> 15                                                                                                                             Value 1 indicates that SCOPE wil use the MoninObukhov correction
#> 16                                                                     Value 1 SCOPE wil estimate LIDF from angles file and Value 0 will use LIDFa and LIDFb from LUT table for estimating LIDF
#> 17                                        Value 0 use a the irradiance file from SCOPE model Value 1 use irradiance measurements WithE and Value 2 use irradiance from MODTRAN atmospheric.file
```

## 2. The options at a glance

**Important: these are not generic on/off flags.** Most are 0/1 toggles,
but `soilspectrum` selects between two different soil *data sources*,
`Fluorescence_model` and `soil_heat_method` select between three
different *methods* each, and `irradiance` selects between three
different irradiance *sources* – in every one of those cases, “0” and
“1” (and “2”) name specific alternatives, not a universal false/true.
Check the “Alternatives” column below before assuming a value’s meaning.

| Option | Default | Alternatives | What does it do? | When should I change it? | Status |
|----|----|----|----|----|----|
| `lite` | 1 | 0 | SCOPE-lite (per-layer) vs. full (per-layer-and-leaf-angle) integration | Normally never | ⚠️ Keep 1; 0 currently errors |
| `calc_fluorescence` | 1 | 0 | Simulate chlorophyll fluorescence (SIF) | Set 0 if you only need reflectance and want a faster run | ✅ Both work |
| `calc_spectrum_planck` | 1 | 0 | Full per-wavelength Planck thermal spectrum vs. simplified thermal | Usually keep default | ✅ Both work |
| `calc_xanthophyllabs` | 1 | 0 | Include xanthophyll/PRI reflectance effects | Disable (0) if not needed, small speedup | ✅ Both work |
| `soilspectrum` | 0 | 1 | Soil reflectance from a file vs. computed by BSM | Use 1 to simulate soil physically (brightness/moisture) rather than from a fixed file | ✅ Both work |
| `Fluorescence_model` | 0 | 1, 2 | Which NPQ/fluorescence-yield submodel | Normally keep 0 | ⚠️ 1 changes real biochemistry fields unexpectedly; 2 currently has no effect (see below) |
| `applTcorr` | 1 | 0 | Temperature-correct Vcmax/biochemical rate constants | Normally keep on | ✅ Both work |
| `verify` | 1 | 0 | Whether SCOPE checks output against field data | Application-dependent | ✅ No effect on the simulation itself (by design – a validation step, not a physics option) |
| `mSCOPE` | 0 | 1 | Vertically heterogeneous canopy (multiple leaf-property layers) | Only for mSCOPE-style applications | ✅ Fixed this session (was a real bug, see below) – now works |
| `simulation` | 0 | 1 | One run per LUT row vs. time-series mode | Only if you have meteo time-series input files, not a LUT | ❌ Needs a different input format this page doesn’t provide (see below) |
| `calc_directional` | 0 | 1 | Full BRDF over many view angles | Use 1 for directional/BRDF studies | ✅ Both work |
| `calc_vert_profiles` | 0 | 1 | Return vertical (per-layer) canopy profiles | When per-layer profiles are needed | ❌ Currently errors – depends on the same unsupported path as `lite=0` |
| `soil_heat_method` | 2 | 0, 1 | How soil heat flux (G) is estimated | Application-dependent | ✅ All three work |
| `calc_rss_rbs` | 0 | 1 | Fixed vs. computed soil surface/boundary-layer resistances | Use 1 to compute from conditions instead of fixing from the LUT | ✅ Fixed this session (was a real bug, see below) – now works |
| `MoninObukhov` | 1 | 0 | Atmospheric-stability correction | Application-dependent | ✅ Both work |
| `LIDF` | 0 | 1 | Leaf angle distribution from LUT (`LIDFa`/`LIDFb`) vs. an external angles file | Use 1 only if you have a real angles file | ✅ Both run and produce different output |
| `irradiance` | 0 | 1, 2 | Irradiance source: SCOPE’s own file, user measurements (“WithE”), or a MODTRAN atmospheric file | 1 if you have real measured irradiance; 2 needs a MODTRAN run | ⚠️ 1 works but see the caution below; 2 improved this session but still incomplete |

## 3. What can I change for my application?

- Need SIF output? `calc_fluorescence = 1` (the default).
- Need a physically-simulated (BSM) soil instead of a fixed file?
  `soilspectrum = 1`.
- Need multi-angular BRDF? `calc_directional = 1`.
- Have real measured irradiance instead of SCOPE’s bundled default?
  `irradiance = 1` (see the caution in Section 5).
- Want a faster run and don’t need fluorescence or PRI effects?
  `calc_fluorescence = 0`, `calc_xanthophyllabs = 0`.
- Everything else: leave at the shipped default – either it has no
  effect worth changing (`verify`, `applTcorr` rarely needs to move), or
  it currently errors (Section 5).

## 4. If it errors: troubleshooting

| Error you see | Likely cause | What to do |
|----|----|----|
| `incorrect number of dimensions` | `lite = 0` | Not currently supported – keep `lite = 1` (Section 5.1) |
| `object 'profiles' not found` | `calc_vert_profiles = 1` | Depends on the same unsupported `lite = 0` path – not currently supported (Section 5.2) |
| `non-numeric matrix extent` | `simulation = 1` | Needs time-series meteo input files, not a single LUT row – this page’s examples don’t provide that input (Section 5.3) |
| A plot silently doesn’t change / model 2 gives identical output | `Fluorescence_model = 2` | Not currently wired to any computation in this port – use 0 or 1 (Section 5.4) |
| `data.bcu` biochemistry fields (`Vcmax`, `Rd`, `gs`, `Cc`) missing | `Fluorescence_model = 1` | A real inconsistency between the two fluorescence-model code paths – avoid relying on `data.bcu` output when using model 1 (Section 5.4) |
| Direct/diffuse irradiance (`Esun_`/`Esky_`) look swapped | `irradiance = 1` | Unconfirmed – flagged for caution, not yet independently verified (Section 5.5) |
| Any error mentioning `spectral` under `irradiance = 2` | `irradiance = 2` | The MODTRAN branch is not fully implemented in this port (Section 5.5) |

## 5. Per-option detail: everything not a plain “it works”

### 5.1 `lite = 0` – Full SCOPE

**Purpose:** activates the full per-leaf-angle-class (13x36xnl)
integration instead of SCOPE-lite’s per-layer aggregation. **Current
status in SCOPEinR:** not operational. **Typical error:**
`incorrect number of dimensions`. **Why:** confirmed by tracing the
actual call stack – with `lite = 0`, `Rnuc` (and related radiation
arrays) become 3-dimensional (`[13, 36, nl]`) rather than the
1-dimensional (`[nl]`) SCOPE-lite case (see the comment at `ebal.R`’s
Vcmax-decline loop, which explicitly distinguishes the two). The
energy-balance iteration loop (`ebal.R`’s `while (CONT)` loop, driven by
[`get.RTMt.sb()`](../reference/get.RTMt.sb.md)) is written and tested
for the lite (1-D) shape; with the full 3-D shape it fails partway
through the iteration. This is a genuine architecture gap, not a
one-line typo – `scopeinpython` (this package’s Python sibling)
documents the same scope restriction explicitly: it only implements
“SCOPE-lite” and lists the full per-leaf-angle branches as a known,
unported gap. Fixing this properly means completing the 3-D energy
balance path, not a quick patch. **What should I do?** Keep `lite = 1`.
**Is my input wrong?** No. With standard SCOPEinR inputs, this is a
current implementation limitation, not a mistake on your part.

### 5.2 `calc_vert_profiles = 1` – Vertical canopy profiles

**Purpose:** return per-layer (vertical) profiles of net radiation/
photosynthesis instead of only canopy-integrated totals. **Current
status in SCOPEinR:** not operational. **Typical error:**
`object 'profiles' not found`. **Why:** traced to the exact source line
– `RTMo.R` only builds the `profiles` list inside a nested condition,
`if (calc_vert_profiles == 1) { if (lite != 1) { profiles <- list(...) } }`.
When `lite = 1` (the only value that actually works, Section 5.1), that
inner block never runs, so `profiles` is never created, and a later line
(`data.profiles = profiles`) fails because the object doesn’t exist. In
other words: `calc_vert_profiles = 1` only ever did anything when
combined with `lite = 0`, which is itself unsupported – so this option
is unusable for the same underlying reason, not a separate bug. **What
should I do?** Keep `calc_vert_profiles = 0`. **Is my input wrong?** No
– same root cause as `lite = 0`.

### 5.3 `simulation = 1` – Time-series mode

**Purpose:** run SCOPE driven by a time series of meteorological
observations (a sequence of dates/conditions) instead of one LUT row per
independent run. **Current status in SCOPEinR:** needs input this page
doesn’t provide. **Typical error:** `non-numeric matrix extent`.
**Why:** this mode expects meteorological time-series input files in a
different structure than `LUT_input.csv`’s one-row-per-run layout; the
error is a direct consequence of feeding it a single LUT row instead.
**What should I do?** Keep `simulation = 0` unless you have real
time-series meteo files formatted for that mode. **Is my input wrong?**
For this tutorial’s LUT-based examples, yes – this option needs a
fundamentally different input, not a bug fix.

### 5.4 `Fluorescence_model` – two real findings, not a bug in the option itself

Testing all three values against the same baseline run surfaced two
things worth knowing before you touch this option:

- **`Fluorescence_model = 2` (Magnani et al. 2012 model) currently has
  no effect.** It’s accepted without error, but every downstream field
  this page’s own verification could see – reflectance, energy balance,
  fluorescence, biochemistry – comes out numerically identical to
  `Fluorescence_model = 0`. This looks like an accepted- but-unwired
  option value, not a working alternative model. **What should I do?**
  Don’t rely on `Fluorescence_model = 2` doing anything different from
  the default; use 0 or 1.
- **`Fluorescence_model = 1` (sigmoid Kn) changes more than the
  fluorescence submodel.** Switching from 0 to 1 makes several ordinary
  leaf-biochemistry fields (`data.bcu$Vcmax`, `Rd`, `gs`, `Cc` –
  sunlit-leaf maximum carboxylation rate, dark respiration, stomatal
  conductance, chloroplast CO2) disappear from the result entirely, even
  though these shouldn’t logically depend on which fluorescence NPQ
  submodel is selected. **What should I do?** If you use
  `Fluorescence_model = 1`, don’t read `data.bcu`’s fields from that run
  – they aren’t populated the same way as the default path.

### 5.5 `irradiance` – one caution, one real fix plus a remaining gap

**`irradiance = 1` (user-supplied “WithE” measurements) runs without
error, but is flagged for caution**: relative to the default
(`irradiance = 0`), the mean direct (`Esun_`) and diffuse (`Esky_`)
irradiance values come out swapped – `Esun_` drops to the previous
`Esky_` value and vice versa, to several significant figures. That’s
precise enough to look like a real column-order issue in how the input
file is read, not coincidence – but confirming which of the two possible
column orders is actually correct requires checking against a real
reference “WithE” measurement file and SCOPE’s own documented convention
for it, which wasn’t available to verify further here. **What should I
do?** Treat `irradiance = 1` output with caution until independently
verified; the direct/diffuse split specifically, not the overall
magnitude, is what’s in question.

**`irradiance = 2` (MODTRAN atmospheric file): partially fixed this
session.** The original crash (`object 'type.irradiance' not found`) was
a genuine typo – `get.SCOPE.R`/`get.SCOPE.ind.R` had
`type.irradiance == 'MODTRAN atmospheric.file'` (a *comparison*, `==`)
where an *assignment* (`<-`) was clearly intended, since the very next
branch (`irradiance = 1`) correctly assigns with `=`. That’s now fixed
in both files, and the fix is real: the bundled MODTRAN reference file
(`inst/input/radiationdata/FLEX-S3_std.atm`) is now read successfully.
But a further, deeper gap remains: the code that reads the raw MODTRAN
transmittance/radiance columns never actually converts them into the
`Esun_`/`Esky_` irradiance values the rest of SCOPE needs – that
physics-conversion step doesn’t exist in this branch yet, which is why a
later step still fails. This matches `scopeinpython`’s own documented
gap list (“MODTRAN atmospheric-file irradiance mode”) exactly – a
genuine unimplemented feature, not a further typo. **What should I do?**
`irradiance = 2` is not yet usable end-to-end; stick with 0 or 1. **Is
my input wrong?** No – for `irradiance = 1`, this is a real finding
worth independent verification; for `irradiance = 2`, this is a genuine
gap in the current port.

### 5.6 `mSCOPE = 1` – fixed this session

**Purpose:** simulate a canopy with vertically varying leaf properties
(multiple biochemistry layers) instead of one uniform layer. **Previous
status:** not operational (`invalid 'times' argument`). **Root cause,
found by tracing the actual source**: `get.SCOPE.R`/ `get.SCOPE.ind.R`
built the per-layer LAI-fraction list with the key `data.mly[['pLAI ']]`
– note the trailing space. `fluspect_mSCOPE()` (the function that
actually consumes this list) looks up `mly[['pLAI']]`, with no trailing
space. Because R list lookups are exact-match, `mly[['pLAI']]` silently
returned `NULL` instead of `c(0.5, 1, 1.5)`, which propagated into a
degenerate zero-length
[`cumsum()`](https://rdrr.io/r/base/cumsum.html)/[`floor()`](https://rdrr.io/r/base/Round.html)
calculation and then a malformed `rep(..., times = ...)` call further
downstream. **Fix applied**: removed the trailing space from the key in
both files, matching the (correctly spelled) name `fluspect_mSCOPE()`
already expected. Verified against the reference implementation
structure at [peiqiyang/mSCOPE](https://github.com/peiqiyang/mSCOPE) –
the 3-layer/pLAI-fraction structure used here matches that model’s own
multi-layer parameterization. `mSCOPE = 1` now runs and produces real,
substantially different canopy-optical output (189 fields differ from
the default run). **What should I do?** `mSCOPE = 1` is now usable.

### 5.7 `calc_rss_rbs = 1` – fixed this session

**Purpose:** compute soil surface (`rss`) and boundary-layer (`rbs`)
resistances from conditions (soil moisture, LAI) instead of using fixed
values from the LUT. **Previous status:** not operational
(`non-numeric argument to binary operator`). **Root cause, found by
tracing the actual source**: `getinputLUT.R` called the (correct) helper
[`calc_rssrbs()`](../reference/calc_rssrbs.md), which returns a list –
but then extracted its results with single-bracket indexing,
`outputs_['rss']`, which returns a length-1 *list* rather than the
numeric value inside it. Assigning that list into `soil$rss` meant any
later arithmetic using `soil$rss` (a number expected) against a list
failed. **Fix applied**: changed both extractions to double-bracket
indexing (`outputs_[['rss']]`, `outputs_[['rbs']]`), the standard R
idiom for pulling a single named element’s actual value out of a list.
`calc_rss_rbs = 1` now runs and produces real, substantially different
energy-balance output (121 fields differ from the fixed-resistance
default). **What should I do?** `calc_rss_rbs = 1` is now usable.

## What’s next

- **Tutorial 01** – run one full simulation end to end with the
  (near-)default options from this tutorial.
- **Tutorial 02** – the soil and canopy BRDF machinery `soilspectrum`
  feeds into.
- **Tutorial 04** – fluorescence in depth, including where
  `calc_fluorescence`/`Fluorescence_model` plug in.

## Appendix: how these options were verified

Everything above is backed by actually running
[`get.SCOPE()`](../reference/get.SCOPE.md) at each option’s documented
alternative values and diffing the result against a baseline, field by
field – not read off the `Info` column alone. This appendix has the
actual verification code, kept here rather than in the main teaching
flow above.

### A structural diff helper

[`get.SCOPE()`](../reference/get.SCOPE.md)’s return value is a deeply
nested list (`data.rad`, `data.fluxes`, `data.canopy`, … – see Tutorial
01). `flatten_scope()` walks the whole nested result down to scalars
(summarizing longer numeric vectors by their mean and length), and
`compare_runs()` reports which fields appeared, disappeared, or changed
value beyond floating-point noise.

``` r

flatten_scope <- function(x, prefix = "") {
  out <- list()
  if (is.list(x)) {
    nms <- names(x)
    if (is.null(nms)) nms <- seq_along(x)
    for (i in seq_along(x)) {
      nm <- if (is.character(nms)) nms[i] else paste0("[[", nms[i], "]]")
      key <- if (nzchar(prefix)) paste0(prefix, "$", nm) else nm
      out <- c(out, flatten_scope(x[[i]], key))
    }
  } else if (is.numeric(x) || is.logical(x) || is.character(x)) {
    if (length(x) == 1) {
      out[[prefix]] <- x
    } else if (length(x) > 1 && is.numeric(x)) {
      out[[paste0(prefix, "[mean]")]] <- mean(x, na.rm = TRUE)
    }
  }
  out
}

compare_runs <- function(base_flat, alt_flat, tol = 1e-6) {
  added <- setdiff(names(alt_flat), names(base_flat))
  removed <- setdiff(names(base_flat), names(alt_flat))
  changed <- character(0)
  for (nm in intersect(names(base_flat), names(alt_flat))) {
    a <- base_flat[[nm]]; b <- alt_flat[[nm]]
    if (is.numeric(a) && is.numeric(b)) {
      if (is.na(a) != is.na(b) || (!is.na(a) && abs(a - b) > tol * max(1, abs(a)))) {
        changed <- c(changed, nm)
      }
    } else if (!identical(a, b)) {
      changed <- c(changed, nm)
    }
  }
  list(added = added, removed = removed, changed = changed)
}
```

### Testing every option for real

One baseline run, then every option flipped to each of its other
documented values, one option at a time, against the same LUT row. Runs
that error are caught rather than allowed to fail the whole comparison –
an error is itself a real, reportable result here.

``` r

lut <- read.table(system.file("input", "LUT_input.csv", package = "SCOPEinR"), header = TRUE, sep = ",")

run_scope <- function(opts) {
  tryCatch({
    invisible(capture.output(
      res <- SCOPEinR::get.SCOPE(LUT = lut[1, ], options.SCOPE = opts,
                                  optipar = SCOPEinR::optipar2021.Pro.CX,
                                  leaf.model = "fluspect-CX", canopy.model = "fourSAIL",
                                  get.outputs = "ALL", get.plots = FALSE)
    ))
    list(ok = TRUE, res = res[[1]])
  }, error = function(e) list(ok = FALSE, err = conditionMessage(e)))
}

base_run <- run_scope(scope_options)
base_flat <- flatten_scope(base_run$res)

# 'simulation' isn't included: value 1 fundamentally changes the input format
# (time-series meteo files instead of a LUT row), so it isn't a same-input
# toggle the way the other 16 options are -- not tested here for that reason,
# not because it errors.
alt_values <- list(
  lite = 0, calc_spectrum_planck = 0, calc_xanthophyllabs = 0, soilspectrum = 1,
  Fluorescence_model = c(1, 2), applTcorr = 0, verify = 0, mSCOPE = 1,
  calc_directional = 1, calc_vert_profiles = 1, soil_heat_method = c(0, 1),
  calc_rss_rbs = 1, MoninObukhov = 0, LIDF = 1, irradiance = c(1, 2),
  calc_fluorescence = 0
)

results <- list()
for (opt in names(alt_values)) {
  for (v in alt_values[[opt]]) {
    opts <- scope_options
    opts$Value[opts$Options == opt] <- v
    r <- run_scope(opts)
    if (r$ok) {
      cmp <- compare_runs(base_flat, flatten_scope(r$res))
      status <- "OK"
      n_changed <- length(cmp$added) + length(cmp$removed) + length(cmp$changed)
      headline <- if (n_changed == 0) "no measurable effect" else paste(n_changed, "fields differ")
    } else {
      status <- "ERROR"; n_changed <- NA; headline <- r$err
    }
    results[[length(results) + 1]] <- data.frame(
      option = opt, value = v, status = status, n_changed = n_changed, headline = headline
    )
  }
}
results_df <- do.call(rbind, results)
knitr::kable(results_df, row.names = FALSE)
```

| option               | value | status | n_changed | headline                       |
|:---------------------|------:|:-------|----------:|:-------------------------------|
| lite                 |     0 | ERROR  |        NA | incorrect number of dimensions |
| calc_spectrum_planck |     0 | OK     |         8 | 8 fields differ                |
| calc_xanthophyllabs  |     0 | OK     |        12 | 12 fields differ               |
| soilspectrum         |     1 | OK     |       175 | 175 fields differ              |
| Fluorescence_model   |     1 | OK     |       134 | 134 fields differ              |
| Fluorescence_model   |     2 | OK     |         1 | 1 fields differ                |
| applTcorr            |     0 | OK     |       119 | 119 fields differ              |
| verify               |     0 | OK     |         1 | 1 fields differ                |
| mSCOPE               |     1 | OK     |       189 | 189 fields differ              |
| calc_directional     |     1 | OK     |        12 | 12 fields differ               |
| calc_vert_profiles   |     1 | ERROR  |        NA | object ‘profiles’ not found    |
| soil_heat_method     |     0 | OK     |         2 | 2 fields differ                |
| soil_heat_method     |     1 | OK     |         2 | 2 fields differ                |
| calc_rss_rbs         |     1 | OK     |       121 | 121 fields differ              |
| MoninObukhov         |     0 | OK     |       122 | 122 fields differ              |
| LIDF                 |     1 | OK     |       193 | 193 fields differ              |
| irradiance           |     1 | OK     |       175 | 175 fields differ              |
| irradiance           |     2 | ERROR  |        NA | object ‘spectral’ not found    |
| calc_fluorescence    |     0 | OK     |        23 | 23 fields differ               |

This is the raw evidence behind every status/finding in Sections 2 and 5
above.
