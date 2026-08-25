# Get Spectral Characteristics for the SPART Model

`get.spectra.spart` Calculates the spectra of hemispherical and
directional observed data. This function retrieves the necessary
spectral characteristics for the SPART model, including spectral ranges
and intervals.

## Usage

``` r
get.spectra.spart(getSpectral = T)
```

## Arguments

- getSpectral:

  A boolean indicating whether to retrieve the spectral object needed
  for SPART (default is TRUE).

## Value

A list containing all spectral ranges and intervals. If `getSpectral` is
TRUE, the list will include the spectral object for SPART.

## Examples

``` r
spectra <- get.spectra.spart(getSpectral = TRUE)
```
