# Get TIFF from a folder and generate the stack.

This function is for Sentinel-2 data using the names B1,B2 ...and SCL

## Usage

``` r
getStacks(rasterFiles = NULL, frequency = "Daily", bands = NULL, output = NULL)
```

## Arguments

- rasterFiles:

  path with the Netcdf

- frequency:

  Daily

- bands:

  a vector with the names of the inputs of the NetCDF

- output:

  path of the outputs

## Value

a stack

## Examples

``` r
if (FALSE) { # \dontrun{
rasterFiles <- "path/to/netcdf/files"
frequency <- "Daily"
bands <- c("B1", "B2", "B3", "SCL")
output <- "path/to/output/folder"
stack <- getStacks(rasterFiles, frequency, bands, output)
} # }
```
