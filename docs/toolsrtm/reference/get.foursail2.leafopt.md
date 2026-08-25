# Leaf optics dispatch used internally by `foursail2`

Computes leaf reflectance/transmittance for a single leaf-parameter row
using whichever of the 5 supported leaf models is requested, mirroring
the dispatch already used in [`foursail`](foursail.md) and
[`inform`](inform.md) (so `foursail2` is no longer limited to
PROSPECT-PRO/PROSPECT-D).

## Usage

``` r
get.foursail2.leafopt(inputRow, LeafModel)
```

## Arguments

- inputRow:

  data.frame. One row of leaf parameters (column names depend on
  `LeafModel`).

- LeafModel:

  character. One of `"PROSPECT-PRO"`, `"PROSPECT-D"`, `"Liberty"`,
  `"Fluspect-B"`, `"Fluspect-B-Cx"`.

## Value

A list whose 2nd element is reflectance and 3rd element is transmittance
(same convention as
`prospect_PRO`/`prospect_DB`/`liberty`/`getFluspect.B`/`getFluspect.Cx`).
