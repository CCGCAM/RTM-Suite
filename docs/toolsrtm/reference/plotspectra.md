# Plot one or more spectra from a data frame

Plot one or more spectra from a data frame

## Usage

``` r
plotspectra(
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

  data.frame. Table containing the spectra to plot, with wavelengths in
  one column and reflectance/radiance values in one or more other
  columns.

- x:

  integer or character. Column index or name in `df` to use as the
  x-axis (typically wavelength). Default 1.

- y:

  integer or character vector. Column index(es) or name(s) in `df` to
  plot as spectra on the y-axis. Default 2.

- cols:

  vector. Line colors, one per element of `y`. Defaults to `y` itself.

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
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) call.

## Value

plot
