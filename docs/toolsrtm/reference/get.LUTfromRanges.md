# Get a LUT based on a table with Min and Max ranges

Get a LUT based on a table with Min and Max ranges

## Usage

``` r
get.LUTfromRanges(
  LUT = NULL,
  nLUT = NULL,
  setseed = 1234,
  leaf.model = "PROSPECT-PRO",
  canopy.model = "fourSAILH",
  distribution = "gauss"
)
```

## Arguments

- LUT:

  a data.frame with three columns (input,min, max)

- nLUT:

  Number of LUT samples to generate

- setseed:

  Random seed for reproducibility

- leaf.model:

  Leaf model to use (e.g., 'PROSPECT-PRO', 'PROSPECT-D', 'Liberty',
  'FLUSPECT-Cx')

- canopy.model:

  Canopy model to use (e.g., 'fourSAILH', 'INFORM')

- distribution:

  Distribution type for LUT generation ('uniform' or 'gauss')

## Value

List containing LUT for given inputs

## Examples

``` r
if (FALSE) { # \dontrun{
# LUT.range: a data.frame with columns (input, min, max) defining the
# sampling range for each PROSPECT/SAIL parameter you want in the LUT —
# build your own with the parameter names your chosen leaf/canopy model
# expects (see the model's own documentation for its parameter names).
# Generate LUT with PROSPECT-PRO and fourSAILH models using a Gaussian distribution
LUT_example <- get.LUTfromRanges(LUT=LUT.range,nLUT = 1000, setseed = 42,
                                leaf.model = 'PROSPECT-PRO',
                                canopy.model = 'fourSAILH',
                                distribution = 'gauss')

# Generate LUT with PROSPECT-D and INFORM models using a Uniform distribution
LUT_example_uniform <- get.LUTfromRanges(LUT=LUT.range,nLUT = 500, setseed = 123,
                                        leaf.model = 'PROSPECT-D',
                                        canopy.model = 'INFORM',
                                        distribution = 'uniform')
} # }

```
