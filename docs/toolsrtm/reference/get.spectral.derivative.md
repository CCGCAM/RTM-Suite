# Calculate spectral derivative

This function calculates spectral derivatives using either finite
approximation or the Savitzky-Golay filter method.

## Usage

``` r
get.spectral.derivative(df, m = 1, method = "sgolay", get.plot = T)
```

## Arguments

- df:

  A DataFrame containing spectral data. Columns represent wavelengths
  and spectral values.

- m:

  The order of the derivative. Default is 1.

- method:

  The method to use for derivative calculation. Options are "finApprox"
  for finite approximation or "sgolay" for Savitzky-Golay filter.
  Default is "sgolay".

- get.plot:

  a boolena TRUE or FALSE; if TRUE a plot will be done

## Value

A DataFrame with spectral derivatives calculated based on the specified
method.

## Examples

``` r
if (FALSE) { # \dontrun{
# Assuming you have a DataFrame named 'spectra_df' with columns representing wavelengths and spectral values
derived_df <- get.spectral.derivative(spectra_df, m = 2, method = "sgolay")
} # }
```
