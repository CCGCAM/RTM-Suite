# latin_hypercube_input function

`latin_hypercube_input` draws a maximin Latin Hypercube sample of
`n_spectra` parameter sets over the parameter ranges given in `tab`, and
writes them to a CSV file (`lh_ts.csv`) for use as SCOPE simulation
input (e.g. to generate a synthetic training/verification dataset). If
`LIDFa` and `LIDFb` are among the sampled variables, the sampled pair is
transformed so that `abs(LIDFa + LIDFb) <= 1` is respected.

## Usage

``` r
latin_hypercube_input(
  tab = read.table(file.path("input", "dataset for_verification", "input_borders.csv"),
    header = TRUE),
  n_spectra = 30,
  outdir = file.path("input", "dataset for_verification")
)
```

## Arguments

- tab:

  data.frame. Parameter bounds table with columns `include` (logical,
  whether the variable is sampled), `lower`/`upper` (parameter bounds),
  and `variable` (parameter name). Defaults to reading
  `input/dataset for_verification/input_borders.csv`.

- n_spectra:

  integer. Number of parameter sets (LHS samples) to draw. Default 30.

- outdir:

  character. Output directory in which `lh_ts.csv` is written (and, if
  `tab` is missing, from which `input_borders.csv` is read). Default
  `input/dataset for_verification`.

## Value

No return value; called for its side effect of writing the sampled
parameter table to `file.path(outdir, 'lh_ts.csv')` and printing a
summary message. Errors if that file already exists.

## Author

    Christiaan van der Tol(Original version in Matlab)

Carlos Camino (Ported version into R)

## Examples

``` r
if (FALSE) { # \dontrun{
latin_hypercube_input(tab, n_spectra = 30, outdir = "input/dataset for_verification")
} # }
```
