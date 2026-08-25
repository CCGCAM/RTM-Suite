# Plot spectra as discrete points at Sentinel band positions

Like [`plotspectra`](plotspectra.md), but plots each spectrum as
discrete points with the x-axis labeled by band name (e.g. Sentinel-2
band names) instead of a continuous wavelength line — for cases where
the input columns are individual sensor bands rather than a continuous
spectrum.

## Usage

``` r
plotspectra_SE(
  df,
  x = 1,
  y = 2,
  cols = y,
  xlim = range(df[, x], na.rm = T),
  ylim = range(df[, y], na.rm = T),
  main = "",
  xlab = "",
  ylab = "",
  ...
)
```

## Arguments

- df:

  data.frame. Table with one column per sensor band plus a values
  column, following the layout produced by
  [`plotspectra_image`](plotspectra_image.md).

- x:

  integer or character. Column index or name in `df` whose column names
  are used as x-axis band labels. Default 1.

- y:

  integer or character vector. Column index(es) or name(s) in `df` to
  plot as points on the y-axis. Default 2.

- cols:

  vector. Point colors, one per element of `y`. Defaults to `y` itself.

- xlim:

  numeric vector of length 2. X-axis plot limits. Defaults to the data
  range of column `x`.

- ylim:

  numeric vector of length 2. Y-axis plot limits. Defaults to the data
  range of columns `y`.

- main:

  character. Plot title.

- xlab:

  character. X-axis label.

- ylab:

  character. Y-axis label.

- ...:

  additional graphical parameters passed to the underlying
  [`points()`](https://rdrr.io/r/graphics/points.html) call.

## Value

a plot
