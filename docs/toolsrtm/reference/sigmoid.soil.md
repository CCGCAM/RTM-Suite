# Calculate Soil Moisture Content Using a Sigmoid Model

This function computes the Soil Moisture Content (SMC) using a sigmoid
model based on the input parameters.

## Usage

``` r
sigmoid.soil(phi, K, a, psi)
```

## Arguments

- phi:

  A numeric value representing the water content parameter.

- K:

  A numeric value representing the sigmoid curve parameter.

- a:

  A numeric value representing the sigmoid curve parameter.

- psi:

  A numeric value representing the sigmoid curve parameter.

## Value

A numeric value representing the Soil Moisture Content (SMC).

## Examples

``` r
phi <- 0.2  # Example water content parameter
K <- 0.5    # Example sigmoid curve parameter
a <- 2.0    # Example sigmoid curve parameter
psi <- 0.1  # Example sigmoid curve parameter
smc <- sigmoid.soil(phi, K, a, psi)  # Calculate Soil Moisture Content
print(smc)  # Output the Soil Moisture Content
#> [1] 0.1688962
```
