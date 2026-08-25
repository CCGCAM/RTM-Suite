# GetIndices: A Function to Estimate Common Spectral Indices

This function calculates common spectral vegetation indices from
reflectance data, using user-specified spectral domains and a scaling
factor if needed.

## Usage

``` r
getIndices(data, pattern.rfl = "R.", factor = NULL, spectral.domain = NULL)
```

## Arguments

- data:

  A dataframe or matrix containing reflectance values in named columns.

- pattern.rfl:

  A character string indicating the prefix pattern for reflectance
  column names. Default is `'R.'`, which assumes column names like
  `'R.400'`, `'R.670'`, etc.

- factor:

  A numeric scaling factor applied to reflectance values (e.g., 1/10000
  for scaling reflectance from digital numbers). Default is `NULL`,
  meaning no scaling is applied (factor = 1).

- spectral.domain:

  A character string specifying the spectral domain for index
  calculation. Available options are:

  - `"VNIR"` (default): Includes wavelengths between 400-850 nm.

  - `"SWIR"`: Includes wavelengths between 800-1750 nm.

  - `"VNIR-SWIR"`: Includes wavelengths between 400-1750 nm. If `NULL`,
    it defaults to `"VNIR"`.

## Value

A dataframe containing the original data along with computed spectral
indices.

## Examples

``` r
# Example 1: Compute spectral indices using default VNIR domain
data_example <- data.frame(R.400 = runif(10, 0, 1), R.670 = runif(10, 0, 1),
                            R.800 = runif(10, 0, 1))
result <- getIndices(data_example)
#> [1] "Estimating indices using VNIR domain..."
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======                                                               |  10%  |                                                                              |==============                                                        |  20%  |                                                                              |=====================                                                 |  30%  |                                                                              |============================                                          |  40%  |                                                                              |===================================                                   |  50%  |                                                                              |==========================================                            |  60%  |                                                                              |=================================================                     |  70%  |                                                                              |========================================================              |  80%  |                                                                              |===============================================================       |  90%  |                                                                              |======================================================================| 100%

# Example 2: Compute indices using SWIR domain with a scaling factor of 1/10000
data_example <- data.frame(R.850 = runif(10, 0, 1), R.1200 = runif(10, 0, 1),
                            R.1600 = runif(10, 0, 1))
result <- getIndices(data_example, spectral.domain = "SWIR", factor = 1/10000)
#> [1] "Estimating indices using SWIR domain..."
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======                                                               |  10%  |                                                                              |==============                                                        |  20%  |                                                                              |=====================                                                 |  30%  |                                                                              |============================                                          |  40%  |                                                                              |===================================                                   |  50%  |                                                                              |==========================================                            |  60%  |                                                                              |=================================================                     |  70%  |                                                                              |========================================================              |  80%  |                                                                              |===============================================================       |  90%  |                                                                              |======================================================================| 100%

# Example 3: Compute indices for VNIR-SWIR domain with custom reflectance column pattern
data_example <- data.frame(Reflectance_450 = runif(10, 0, 1), Reflectance_900 = runif(10, 0, 1))
result <- getIndices(data_example, pattern.rfl = "Reflectance_", spectral.domain = "VNIR-SWIR")
#> [1] "Estimating indices using VNIR-SWIR domain..."
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======                                                               |  10%  |                                                                              |==============                                                        |  20%  |                                                                              |=====================                                                 |  30%  |                                                                              |============================                                          |  40%  |                                                                              |===================================                                   |  50%  |                                                                              |==========================================                            |  60%  |                                                                              |=================================================                     |  70%  |                                                                              |========================================================              |  80%  |                                                                              |===============================================================       |  90%  |                                                                              |======================================================================| 100%
```
