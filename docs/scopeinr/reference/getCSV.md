# Get CSV Data from Multiple Folders

This function reads CSV files from multiple folders and combines them
into a list.

## Usage

``` r
getCSV(path.out = NULL, n.folders = 4, files.names = "All")
```

## Arguments

- path.out:

  The path to the root folder containing the subfolders with CSV files.

- n.folders:

  The number of folders from the end of the directory structure to
  consider.

- files.names:

  The type of files to look for. Options: "Fluorescence", "Reflectance",
  or "Radiance".

## Value

A list with `data` (list, one element per file pattern in `files.names`,
each the row-bound data across the last `n.folders` folders), `names`
(the file patterns with the `.csv` suffix stripped), and `lut` (the
row-bound `Parameters/inputLUT.csv` across those same folders).

## Examples

``` r
if (FALSE) { # \dontrun{
path.out <- "outs"
n.folders <- 4
files.names <- "Fluorescence"
combined_data <- getCSV(path.out, n.folders, files.names)
} # }
```
