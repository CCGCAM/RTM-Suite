# Cost function for Fluspect parameter fitting

Cost function for Fluspect parameter fitting

## Usage

``` r
COST_4Fluspect(params, measurement, input)
```

## Arguments

- params:

  numeric vector of length 8, in order: Cab, Cdm, Cw, Cs, Cca, Cant, Cx,
  N - leaf biochemistry/structure parameters being fitted.

- measurement:

  the observed (measured) spectrum being fitted against.

- input:

  list of 6 elements, in order: `leafbio` (baseline leaf parameters),
  `optipar` (optical parameters), `spectral` (spectral configuration),
  `include` (0/1 flags per parameter, controlling whether each is fitted
  or held at its `leafbio` default), `target`, `range`.

## Value

rfl and tran resampled
