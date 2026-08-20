# Calculate Soil Moisture Content Using a Sigmoid Model

This function computes the Soil Moisture Content (SMC) using a sigmoid
model based on the input parameters.

This function calculates a sigmoid curve based on input parameters.

## Usage

``` r
sigmoid.soil(x, K, a, psi)

sigmoid.soil(x, K, a, psi)
```

## Arguments

- x:

  Numeric. Input value.

- K:

  Numeric. Maximum value of the sigmoid.

- a:

  Numeric. Steepness of the curve.

- psi:

  Numeric. Offset parameter.

- phi:

  A numeric value representing the water content parameter.

## Value

A numeric value representing the Soil Moisture Content (SMC).

Numeric. Sigmoid value for the input x.

## Examples

``` r
phi <- 0.2  # Example water content parameter
K <- 0.5    # Example sigmoid curve parameter
a <- 2.0    # Example sigmoid curve parameter
psi <- 0.1  # Example sigmoid curve parameter
smc <- sigmoid.soil(phi, K, a, psi)  # Calculate Soil Moisture Content
print(smc)  # Output the Soil Moisture Content
#> [1] 0.1688962

if (FALSE) { # \dontrun{
# Example of sigmoid function
result <- sigmoid.soil(0.5, K = 1, a = 5, psi = 2)
print(result)

} # }
```
