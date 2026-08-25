
# Define the get_spectrum function for MARMIT-1
#' Calculate Spectral Reflectance for Wet Soils Using MARMIT-1
#'
#' This function computes the spectral reflectance of wet soils using the MARMIT-1 model.
#'
#' @param n A numeric vector of the spectral optical index of water.
#' @param alpha A numeric vector of water absorption spectral coefficient in cm^-1.
#' @param Rd A numeric vector of the reflectance of dry soil.
#' @param L A numeric value representing the thickness of the water layer in cm.
#' @param eps A numeric value (between 0 and 1) representing the wet soil surface ratio.
#'
#' @return A numeric vector representing the reflectance of wet soil.
#' @export
#'
#' @examples
#' n <- seq(1.33, 1.40, length.out = 100)  # Example optical index of water
#' alpha <- seq(0.1, 0.2, length.out = 100)  # Example water absorption coefficient
#' Rd <- rep(0.3, 100)  # Example reflectance of dry soil
#' L <- 0.01  # Water layer thickness in cm
#' eps <- 0.5  # Wet soil surface ratio
#' Rm <- get.marmit1(n, alpha, Rd, L, eps)  # Calculate wet soil reflectance
#' plot(Rm, type = "l", col = "blue", xlab = "Wavelength Index", ylab = "Reflectance")
#'
get.marmit1 <- function(n, alpha, Rd, L, eps) {

  # Calculate r12_diffuse
  r12_diffuse <- ((3 * n^2 + 2 * n + 1) / (3 * (n + 1)^2)
                  - 2 * n^3 * (n^2 + 2 * n - 1) / ((n^2 + 1)^2 * (n^2 - 1))
                  + n^2 * (n^2 + 1) * log(n) / (n^2 - 1)^2
                  - n^2 * (n^2 - 1)^2 * log(n * (n + 1) / (n - 1)) / (n^2 + 1)^3)

  # Calculate t12_diffuse and t21_diffuse
  t12_diffuse <- 1 - r12_diffuse
  r21_diffuse <- 1 - (1 - r12_diffuse) / n^2
  t21_diffuse <- 1 - r21_diffuse

  # Calculate Tw_diffuse
  if (L > 0) {
    Tw_diffuse <- (1 - alpha * L) * exp(-alpha * L) + (alpha * L)^2 * exp1_base(alpha * L)
  } else {
    Tw_diffuse <- rep(1, length(n))  # If L = 0, set Tw_diffuse to 1
  }

  # MARMIT model with diffuse light in water layer
  Rw <- (t12_diffuse * t21_diffuse * Rd * Tw_diffuse^2) / (1 - r21_diffuse * Rd * Tw_diffuse^2)
  Rm <- eps * Rw + (1 - eps) * Rd  # Wet soil reflectance

  return(Rm)
}
# Define the sigmoid function
#' Calculate Soil Moisture Content Using a Sigmoid Model
#'
#' This function computes the Soil Moisture Content (SMC) using a sigmoid model based on the input parameters.
#'
#' @param phi A numeric value representing the water content parameter.
#' @param K A numeric value representing the sigmoid curve parameter.
#' @param a A numeric value representing the sigmoid curve parameter.
#' @param psi A numeric value representing the sigmoid curve parameter.
#'
#' @return A numeric value representing the Soil Moisture Content (SMC).
#' @export
#'
#' @examples
#' phi <- 0.2  # Example water content parameter
#' K <- 0.5    # Example sigmoid curve parameter
#' a <- 2.0    # Example sigmoid curve parameter
#' psi <- 0.1  # Example sigmoid curve parameter
#' smc <- sigmoid.soil(phi, K, a, psi)  # Calculate Soil Moisture Content
#' print(smc)  # Output the Soil Moisture Content
#'
sigmoid.soil <- function(phi, K, a, psi) {

  SMC <- K / (1 + a * exp(-psi * phi))  # Soil Moisture Content
  return(SMC)
}

# Define a function to compute the exponential integral E1(x)
#' Compute Exponential Integral E1(x)
#'
#' This function computes the exponential integral E1(x) for each element in the input vector.
#'
#' @param x A numeric vector of values.
#'
#' @return A numeric vector containing the result of the exponential integral E1(x) for each value in `x`.
#' @export
#'
#' @examples
#' exp1_base(1:10)  # Compute E1 for values 1 to 10
#'
exp1_base <- function(x) {
  sapply(x, function(xi) {
    integrate(function(t) exp(-t) / t, lower = xi, upper = Inf)$value
  })
}
