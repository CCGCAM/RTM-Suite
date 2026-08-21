#' Simulation of Spectral Reflectance of Wet Soils in the Solar Domain using MARMIT-2
#'
#' This function calculates the reflectance of wet soils using the MARMIT-2 model
#' based on various input parameters such as the optical index of water, refractive index
#' of soil particles, and dry soil reflectance.
#'
#' @param n_w Numeric vector. Spectral optical index of water.
#' @param alpha_w Numeric vector. Spectral specific absorption of water in cm-1.
#' @param n_i Numeric. Real part of the refractive index of soil particles.
#' @param k_i Numeric. Imaginary part of the refractive index of soil particles.
#' @param Rd Numeric vector. Reflectance of dry soil.
#' @param L Numeric. Thickness of the water layer.
#' @param eps Numeric. Fraction of the surface that is wet.
#' @param d_i Numeric. Fraction of soil in the mixture.
#' @param wls Numeric vector. Wavelengths in nanometers (nm).
#'
#' @return Numeric vector. Reflectance of wet soil.
#'
#' @examples
#' \dontrun{
#' # Example input data
#' n_w <- c(1.33, 1.34, 1.35)
#' alpha_w <- c(0.01, 0.015, 0.02)
#' n_i <- 1.5
#' k_i <- 0.02
#' Rd <- c(0.3, 0.35, 0.4)
#' L <- 0.1
#' eps <- 0.6
#' d_i <- 0.5
#' wls <- c(400, 500, 600)
#'
#' # Calculate wet soil reflectance
#' Rm <- get_spectrum(n_w, alpha_w, n_i, k_i, Rd, L, eps, d_i, wls)
#' print(Rm)
#'
#' }
get.marmit2 <- function(n_w, alpha_w, n_i, k_i, Rd, L, eps, d_i, wls) {


  # Imaginary part of refractive index of water
  k_w <- alpha_w * wls * 10^-7 / (4 * pi)
  # Complex permittivity of water
  e_w <- (n_w + 1i * k_w)^2
  # Complex permittivity of soil particles
  e_i <- (n_i + 1i * k_i)^2
  # Dielectric average
  e <- d_i * e_i + (1 - d_i) * e_w
  # Effective refractive index of the mixture
  n <- Re(sqrt(e))
  k <- Im(sqrt(e))
  # Effective absorption coefficient
  alpha <- 4 * pi * k / (wls * 10^-7)
  # Fresnel coefficients integrated over the hemisphere
  r12_diffuse <- (3 * n^2 + 2 * n + 1) / (3 * (n + 1)^2) -
    2 * n^3 * (n^2 + 2 * n - 1) / ((n^2 + 1)^2 * (n^2 - 1)) +
    n^2 * (n^2 + 1) * log(n) / (n^2 - 1)^2 -
    n^2 * (n^2 - 1)^2 * log(n * (n + 1) / (n - 1)) / (n^2 + 1)^3

  t12_diffuse <- 1 - r12_diffuse
  r21_diffuse <- 1 - (1 - r12_diffuse) / n^2
  t21_diffuse <- 1 - r21_diffuse
  # Transmission of the water layer
  if (L > 0) {
    Tw_diffuse <- (1 - alpha * L) * exp(-alpha * L) + (alpha * L)^2 *  exp1_base(alpha * L)
  } else {
    Tw_diffuse <- rep(1, length(n))
  }
  # Reflectance of totally wet soil
  Rw <- t12_diffuse * t21_diffuse * Rd * Tw_diffuse^2 / (1 - r21_diffuse * Rd * Tw_diffuse^2)
  # Mixing the reflectances of wet and dry areas
  Rm <- (eps * Rw^(1 / 2.27) + (1 - eps) * Rd^(1 / 2.27))^2.27

  return(Rm)
}

# sigmoid.soil() is defined once, in marmit1.R (used here via get.marmit.rsoil()).

# Define a function to compute the exponential integral E1(x)
exp1_base <- function(x) {
  sapply(x, function(xi) {
    integrate(function(t) exp(-t) / t, lower = xi, upper = Inf)$value
  })
}
