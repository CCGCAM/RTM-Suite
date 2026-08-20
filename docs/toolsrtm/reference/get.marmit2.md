# Simulation of Spectral Reflectance of Wet Soils in the Solar Domain using MARMIT-2

This function calculates the reflectance of wet soils using the MARMIT-2
model based on various input parameters such as the optical index of
water, refractive index of soil particles, and dry soil reflectance.

## Usage

``` r
get.marmit2(n_w, alpha_w, n_i, k_i, Rd, L, eps, d_i, wls)
```

## Arguments

- n_w:

  Numeric vector. Spectral optical index of water.

- alpha_w:

  Numeric vector. Spectral specific absorption of water in cm-1.

- n_i:

  Numeric. Real part of the refractive index of soil particles.

- k_i:

  Numeric. Imaginary part of the refractive index of soil particles.

- Rd:

  Numeric vector. Reflectance of dry soil.

- L:

  Numeric. Thickness of the water layer.

- eps:

  Numeric. Fraction of the surface that is wet.

- d_i:

  Numeric. Fraction of soil in the mixture.

- wls:

  Numeric vector. Wavelengths in nanometers (nm).

## Value

Numeric vector. Reflectance of wet soil.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example input data
n_w <- c(1.33, 1.34, 1.35)
alpha_w <- c(0.01, 0.015, 0.02)
n_i <- 1.5
k_i <- 0.02
Rd <- c(0.3, 0.35, 0.4)
L <- 0.1
eps <- 0.6
d_i <- 0.5
wls <- c(400, 500, 600)

# Calculate wet soil reflectance
Rm <- get_spectrum(n_w, alpha_w, n_i, k_i, Rd, L, eps, d_i, wls)
print(Rm)

} # }
```
