
#' Define constant for FLuxPesp
#'
#' @return constants
#' @export
#'
define_constants <- function() {
  const <- list(
    A = 6.02214E23, # [mol-1]       Constant of Avogadro
    h = 6.6262E-34, # [J s]         Planck's constant
    c = 299792458,  # [m s-1]       Speed of light
    cp = 1004,      # [J kg-1 K-1]  Specific heat of dry air
    R = 8.314,      # [J mol-1K-1]  Molar gas constant
    rhoa = 1.2047,  # [kg m-3]      Specific mass of air
    g = 9.81,       # [m s-2]       Gravity acceleration
    kappa = 0.4,    # []            Von Karman constant
    MH2O = 18,      # [g mol-1]     Molecular mass of water
    Mair = 28.96,   # [g mol-1]     Molecular mass of dry air
    MCO2 = 44,      # [g mol-1]     Molecular mass of carbon dioxide
    sigmaSB = 5.67E-8, # [W m-2 K-4] Stefan Boltzman constant
    deg2rad = pi/180,  # [rad]       Conversion from deg to rad
    C2K = 273.15       # [K]         Melting point of water
  )
  return(const)
}


#' define spectral region
#'
#' @return spectral regions
#' @export
#'
define_bands <- function() {
  spectral <- list()
  
  # Define spectral regions for SCOPE v_1.40
  # All spectral regions are defined here as row vectors
  # WV Jan. 2013
  
  # 3 spectral regions for SCOPE
  
  reg1 <-  seq(400,2400,1)
  reg2 <- seq(2500,15000,100)
  reg3 <- seq(16000,50000,1000)
  
  spectral$wlS  <- c(reg1, reg2, reg3)
  
  # Other spectral (sub)regions
  spectral$wlP   <- reg1 # PROSPECT data range
  spectral$wlE   <- as.integer(seq(400, 750, 1))
  spectral$wlF   <- as.integer(seq(640, 850, 1)) # chlorophyll fluorescence in E-F matrix
  spectral$wlO   <- reg1 # optical part
  spectral$wlT   <- c(reg2, reg3) # thermal part
  
  # Other spectral (sub)regions
  
  spectral$wlP   <- reg1                            # PROSPECT data range
  spectral$wlE   <- seq(400,750,1)                       # excitation in E-F matrix
  spectral$wlF   <- seq(640,850,1)                       # chlorophyll fluorescence in E-F matrix
  spectral$wlO   <- reg1                            # optical part
  spectral$wlT   <- c(reg2, reg3)                   # thermal part
  
  wlS            <- spectral$wlS
  spectral$wlPAR <- wlS[wlS >= 400 & wlS <= 700]     # PAR range
  
  # Data used by aggreg routine to read in MODTRAN data
  
  spectral$SCOPEspec$nreg  <- 3
  spectral$SCOPEspec$start <- c(400, 2500, 16000)
  spectral$SCOPEspec$end   <- c(2400, 15000, 50000)
  spectral$SCOPEspec$res   <- c(1, 100, 1000)
  
  return(spectral)
}



#' @title num-jacobian function
#'
#' @param x numeric vector. Leaf biochemistry/structure parameters at which to compute the numerical Jacobian.
#' @param spectral list. Spectral configuration (wavelength ranges/resolution) used by the Fluspect model.
#' @param inputLeaf list. Leaf input parameters passed through to \code{calc_fluspect_bcar}.
#' @param optipar list. Optical parameters (refractive index, specific absorption coefficients) used by the Fluspect model.
#'
#' @author Wout Verhoef, Christiaan van der Tol, Joris Timmermans, 
#' @author Ported to R: Carlos. Camino
#' 
#' @return numjacobian  values
#' @export
numjacobian <- function(x, spectral, inputLeaf, optipar) {
  n <- length(x)
  res <- ToolsRTM::calc_fluspect_bcar(x, spectral, leafbio, optipar)
  fx <- c(res$refl, res$tran)
  step <- 1e-6
  J <- array(0, dim = c(length(fx), 2, n))
  
  for (i in 1:n) {
    xstep <- x
    xstep[i] <- x[i] + step
    res_step <- ToolsRTM::calc_fluspect_bcar(xstep, spectral, inputLeaf, optipar)
    fxstep <- c(res_step$refl, res_step$tran)
    J[, , i] <- (fxstep - fx) / step
  }
  
  return(J)
}

#' @title calc_fluspect_bcar function
#'
#' @param params numeric vector of length 6, in order: Cab (chlorophyll content), Cdm (dry matter content), Cw (water content), Cs (senescent material/brown pigments), Cca (carotenoid content), N (leaf structure parameter).
#' @param spectral list. Spectral configuration (wavelength ranges/resolution) used by the Fluspect model.
#' @param leafbio list. Leaf biochemistry/structure parameters; \code{Cab}, \code{Cdm}, \code{Cw}, \code{Cs}, \code{Cca}, and \code{N} are overwritten from \code{params} before simulation.
#' @param optipar list. Optical parameters (refractive index, specific absorption coefficients) used by the Fluspect model.
#'
#' @return rfl and trans
#' @export
#'
calc_fluspect_bcar <- function(params, spectral, leafbio, optipar) {
  # The rest of this function's body (below the stop()) called a function,
  # fluspect_B_CX_PSI_PSII_combined(), that does not exist anywhere in this
  # package -- removed as unreachable dead code (R CMD check's static
  # analysis flags a reference to an undefined global even in unreachable
  # code). See the stop() message for why this can't just be rewired to
  # getFluspect.Cx()/getFluspect.B(): different, unconfirmed parameterization.
  stop(
    "calc_fluspect_bcar() is currently broken: it called fluspect_B_CX_PSI_PSII_combined(), ",
    "which does not exist in this package. The existing getFluspect.Cx()/getFluspect.B() use a ",
    "different parameter set (EWT, LMA, Cx...) than this function's leafbio (Cdm, Cw, Cca...), ",
    "so they can't be substituted without confirming the correct physical mapping between the two ",
    "parameterizations. Needs review before use - see ToolsRTM issue tracker / ask the package author."
  )
}


#' Function for calculating the energy of a single photon given its wavelength.
#'
#' @param lambda numeric. Photon wavelength, in meters.
#'
#' @return e
#' @export
#'
ephoton <- function(lambda) {
  h <- 6.62607004E-34 # Planck constant [J.s]
  c <- 2.99792458E8 # Speed of light [m/s]
  e <- h * c / lambda
  return(e)
}

#' @title  Simpson integration
#'
#' @param y  must be any vector (rows, columns), but of the same length
#' @param x  must be any vector (rows, columns), but of the same length; must be a monotonically increasing series
#' develope by WV Jan. 2013, for SCOPE 1.40
#' @return int
#' @export
#'
#' 
Sint <- function(y, x) {

  nx <- length(x)
  if (length(dim(x)) == 1) {
    x <- t(x)
  }
  if (length(dim(y)) != 1) {
    y <- t(y)
  }
  step <- x[2:nx] - x[1:(nx-1)]
  mean <- 0.5 * (y[1:(nx-1)] + y[2:nx])
  int <- mean * step
  return(int)
}



#' @title Cost function for Fluspect parameter fitting
#'
#' @param params numeric vector of length 8, in order: Cab, Cdm, Cw, Cs, Cca, Cant, Cx, N - leaf biochemistry/structure parameters being fitted.
#' @param measurement the observed (measured) spectrum being fitted against.
#' @param input list of 6 elements, in order: \code{leafbio} (baseline leaf parameters), \code{optipar} (optical parameters), \code{spectral} (spectral configuration), \code{include} (0/1 flags per parameter, controlling whether each is fitted or held at its \code{leafbio} default), \code{target}, \code{range}.
#'
#' @return rfl and tran resampled
#' @export
#'
COST_4Fluspect <- function(params, measurement, input) {
  # Rest of the body removed as unreachable dead code -- same reason as
  # calc_fluspect_bcar() above, which this function calls.
  stop(
    "COST_4Fluspect() is currently broken: it depends on calc_fluspect_bcar(), which called ",
    "fluspect_B_CX_PSI_PSII_combined() - a function that does not exist in this package. ",
    "Needs review before use - see ToolsRTM issue tracker / ask the package author."
  )
}

