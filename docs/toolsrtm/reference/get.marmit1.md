# Calculate Spectral Reflectance for Wet Soils Using MARMIT-1

This function computes the spectral reflectance of wet soils using the
MARMIT-1 model.

## Usage

``` r
get.marmit1(n, alpha, Rd, L, eps)
```

## Arguments

- n:

  A numeric vector of the spectral optical index of water.

- alpha:

  A numeric vector of water absorption spectral coefficient in cm^-1.

- Rd:

  A numeric vector of the reflectance of dry soil.

- L:

  A numeric value representing the thickness of the water layer in cm.

- eps:

  A numeric value (between 0 and 1) representing the wet soil surface
  ratio.

## Value

A numeric vector representing the reflectance of wet soil.

## Examples

``` r
n <- seq(1.33, 1.40, length.out = 100)  # Example optical index of water
alpha <- seq(0.1, 0.2, length.out = 100)  # Example water absorption coefficient
Rd <- rep(0.3, 100)  # Example reflectance of dry soil
L <- 0.01  # Water layer thickness in cm
eps <- 0.5  # Wet soil surface ratio
Rm <- get.marmit1(n, alpha, Rd, L, eps)  # Calculate wet soil reflectance
plot(Rm, type = "l", col = "blue", xlab = "Wavelength Index", ylab = "Reflectance")

```
