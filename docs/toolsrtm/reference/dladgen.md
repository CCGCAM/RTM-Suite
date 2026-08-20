# Computes the leaf angle distribution function value (freq)

### Using the original bimodal distribution function initially proposed in SAIL References

(Verhoef1998) Verhoef, Wout. Theory of radiative transfer models applied
in optical remote sensing of vegetation canopies. Nationaal Lucht en
Ruimtevaartlaboratorium, 1998.
http://library.wur.nl/WebQuery/clc/945481.

## Usage

``` r
dladgen(a, b)
```

## Arguments

- a:

  controls the average leaf slope

- b:

  controls the distribution's bimodality LIDF type a b Planophile 1 0
  Erectophile -1 0 Plagiophile 0 -1 Extremophile 0 1 Spherical -0.35
  -0.15 Uniform 0 0 requirement: \|LIDFa\| + \|LIDFb\| \< 1

## Value

foliar_distrib list. lidf and litab
