# Field data expected here

`1-GetSimulationsLUTs.R` (step 5 onward), `1-ReduceLUTs.R`, `2-Model.R`, and
`2-Model_fast.R` compare/calibrate PROSAIL simulations against a real field
campaign: Carlos Camino's 2019/2020 Vcmax field dataset. That file is private
research data and is **not bundled with this repo**.

To run these scripts, place your own copy here as:

```
2019_2020RawData_Vcmax.csv
```

Expected columns:
- Hyperspectral reflectance columns named `X<wavelength>` (e.g. `X450`, `X451`, ...)
- A `Vcmax` column (measured maximum carboxylation rate)
- Any other field metadata columns used downstream (e.g. `Plot`)

`1-GetSimulationsLUTs.R` resamples this to sensor bands and writes
`2019_2020RawData_Vcmax_withRFL_resampled.csv` to `outs/ForPROSAIL/FieldData/`,
which is what the other three scripts read as their own field-data input.

Without this file, only steps 1-4 of `1-GetSimulationsLUTs.R` (LUT generation,
soil mixing, PROSAIL simulation, sensor resampling) can run — everything that
doesn't need a field comparison.
