# Run SCOPE simulations in parallell

Run SCOPE simulations in parallell

## Usage

``` r
get.SCOPE.parallel(
  LUT,
  options.SCOPE = NULL,
  optipar = NULL,
  path.out = NULL,
  parallel = TRUE,
  canopy.model = "fourSAIL",
  leaf.model = "fluspect-CX",
  get.outputs = "Main",
  get.plots = TRUE,
  get.csv = T,
  n.cores = 8
)
```

## Arguments

- LUT:

  Leaf, biochemistry, viewing angles, meteo, and canopy properties
  needed for running SCOPE model.

- options.SCOPE:

  Optical leaf properties, the total irradiance if a specific value is
  provided instead of the usual Modtran output.

- optipar:

  list. Leaf-level optical parameters passed on to
  [`get.SCOPE`](get.SCOPE.md) for each LUT chunk (e.g.
  [`SCOPEinR::optipar2021.Pro.CX`](optipar-datasets.md)).

- path.out:

  Folder for saving the SCOPE outputs.

- parallel:

  Logical, indicating whether to use parallel processing. Default is
  TRUE.

- canopy.model:

  Selection of canopy model options available are 'fourSAIL', 'fourSAIL'
  and 'INFORM'. By default, the fourSAIL model will be used.

- leaf.model:

  Selection of canopy model options available are 'fluspect-CX',
  'fluspect-B', 'PROSPECT', 'Liberty'. By default, the fluspect-CX model
  will be used.

- get.outputs:

  If get.outputs = 'ALL', all variables will be retrieved; if
  get.outputs = 'Main', only the main variables will be retrieved. By
  default, SCOPE uses get.outputs = 'ALL', for processing a huge LUT, it
  is recommended to use get.outputs = 'Main'.

- get.plots:

  Logical, indicating whether to plot the intermediate plots. Default is
  TRUE.

- get.csv:

  logical, indicating whether to save outputs. Default is TRUE.

- n.cores:

  Integer, indicating the number of cores to use. if this parameter is
  null or missing Detect cores - 2 will be used

## Value

A list of SCOPE simulations.

## Examples

``` r
if (FALSE) { # \dontrun{
get.SCOPE.parallel(LUT, options.SCOPE = NULL, optipar = NULL,
                   path.out = "output", parallel = TRUE, leaf.model = 'fluspect-CX',
                   canopy.model = 'fourSAIL', get.outputs = 'Main', get.plots = TRUE)
} # }
```
