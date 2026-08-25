# Build a LUT with a per-trait distribution choice (Uniform or Gaussian) and an optional Car~Cab correlation

Alternative to [`getLUT`](getLUT.md) for when different traits need
different sampling distributions in the same LUT (e.g. LAI sampled
Uniform while Cab is sampled Gaussian), rather than one distribution for
every trait. Gaussian traits are drawn by truncated-normal rejection
sampling ([`gauss_byMin_Max`](gauss_byMin_Max.md)), so they still
respect `minval`/`maxval` bounds.

## Usage

``` r
get_distributionLUT(
  minval = NULL,
  maxval = NULL,
  nSamples = NULL,
  TypeDistrib = NULL,
  Mean_gauss = NULL,
  Std_gauss = NULL,
  DepCab = NULL,
  setseed = NULL
)
```

## Arguments

- minval:

  one-row data.frame/list, column names = trait names, minimum value per
  trait.

- maxval:

  one-row data.frame/list, column names = trait names, maximum value per
  trait.

- nSamples:

  numeric. Number of LUT rows to generate.

- TypeDistrib:

  named list, one entry per trait in `minval`, each either `"Uniform"`
  or `"Gaussian"`.

- Mean_gauss:

  one-row data.frame/list, mean per trait – only read for traits where
  `TypeDistrib` is `"Gaussian"`.

- Std_gauss:

  one-row data.frame/list, standard deviation per trait – only read for
  traits where `TypeDistrib` is `"Gaussian"`.

- DepCab:

  logical. If `TRUE` and `"Car"` is one of the traits, Car is not drawn
  independently – it's redrawn as
  [`correlatedValue`](correlatedValue.md)`(Cab/4, r = 0.8)`, the
  empirical Cab-Car co-variation seen in leaf pigment data.

- setseed:

  integer. Random seed.

## Value

LUT as a data.frame, one column per trait in `minval`.

## Examples

``` r
minv <- data.frame(Cab = 10, Car = 2, LAI = 0.5)
maxv <- data.frame(Cab = 80, Car = 20, LAI = 7)
distrib <- list(Cab = "Gaussian", Car = "Uniform", LAI = "Uniform")
meang <- data.frame(Cab = 40, Car = NA, LAI = NA)
stdg  <- data.frame(Cab = 15, Car = NA, LAI = NA)
LUT <- get_distributionLUT(minval = minv, maxval = maxv, nSamples = 100,
                            TypeDistrib = distrib, Mean_gauss = meang, Std_gauss = stdg,
                            DepCab = TRUE, setseed = 1)
```
