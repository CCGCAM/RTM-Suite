# Generate LUT for SCOPE (v1)

This function generates samples of a Look-Up Table (LUT) to analyze
variations in plant traits required for a specific scope.

## Usage

``` r
getLUT.SCOPE.v1(inputLUT, nLUT = 100, setseed = 123)
```

## Arguments

- inputLUT:

  A table (LUT) containing data.

- nLUT:

  Number of samples to generate from the input LUT. Defaults to 100.

- setseed:

  Seed value to ensure reproducibility of results. Defaults to 123.

## Value

Returns a data frame containing the generated LUT samples.

## Examples

``` r
if (FALSE) { # \dontrun{
inputLUT <- read.table(system.file("input", "inputs_SCOPE.csv", package = "SCOPEinR"),
                       header = TRUE, sep = ",")
getLUT.SCOPE.v1(inputLUT = inputLUT, nLUT = 100, setseed = 123)
} # }
```
