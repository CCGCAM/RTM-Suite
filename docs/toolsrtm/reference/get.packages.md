# Get the packages

This function retrieves the list of R package dependencies for the
specified path using
[`renv::dependencies()`](https://rstudio.github.io/renv/reference/dependencies.html).
If `renv` is not installed, it installs `renv` first. THis is designed
for getting the packages that the Shiny app needs in this package

## Usage

``` r
get.packages(path = ".")
```

## Arguments

- path:

  The directory path to check for dependencies. Defaults to the current
  directory.

## Value

A character vector of unique package names found in the project
dependencies.

## Examples

``` r
if (FALSE) { # \dontrun{
# Get the packages for the current directory
get.packages()

# Get the packages for a specified path
get.packages("path/to/project")

} # }
```
