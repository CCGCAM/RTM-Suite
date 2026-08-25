
#' Define physical constants used by SCOPEinR
#'
#' \code{define.constants} returns a list of the physical constants used
#' throughout the SCOPEinR model (radiative transfer, energy balance and
#' biochemical sub-models), analogous to the \code{SCOPEinR::constants} table.
#'
#' @return A list with named numeric elements: \code{A} (Avogadro's number, mol^-1), \code{h} (Planck's constant, J s), \code{c} (speed of light, m s^-1), \code{cp} (specific heat of dry air, J kg^-1 K^-1), \code{R} (molar gas constant, J mol^-1 K^-1), \code{rhoa} (specific mass of air, kg m^-3), \code{g} (gravitational acceleration, m s^-2), \code{kappa} (Von Karman constant, dimensionless), \code{MH2O} (molecular mass of water, g mol^-1), \code{Mair} (molecular mass of dry air, g mol^-1), \code{MCO2} (molecular mass of CO2, g mol^-1), \code{sigmaSB} (Stefan-Boltzmann constant, W m^-2 K^-4), \code{deg2rad} (degrees-to-radians conversion factor), \code{C2K} (Celsius-to-Kelvin offset, K).
#' @export
#'
#' @author 	Wout Verhoef (Original version in Matlab)
#' @author Carlos Camino (Ported version into R)
#'
#' @examples
#' const <- define.constants()
define.constants <- function() {
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


#' Define spectral regions used by SCOPEinR
#'
#' \code{define.bands} defines the wavelength grid and spectral sub-regions
#' used throughout SCOPEinR (optical, thermal, PAR, excitation-fluorescence
#' matrix, etc.), following the SCOPE v1.40 spectral band definition (three
#' regions of increasing resolution: 400-2400 nm at 1 nm, 2500-15000 nm at
#' 100 nm, and 16000-50000 nm at 1000 nm).
#'
#' @return A list \code{spectral} with (among others): \code{wlS} (full wavelength vector, nm), \code{wlP} (PROSPECT range), \code{wlE} (excitation wavelengths for the E-F fluorescence matrix), \code{wlF} (chlorophyll fluorescence emission wavelengths), \code{wlO} (optical part), \code{wlT} (thermal part), \code{wlPAR} (PAR range, 400-700 nm), and \code{SCOPEspec} (region boundaries and resolutions used by \code{\link{aggreg}} to read MODTRAN data).
#' @export
#'
#' @author 	Wout Verhoef (Original version in Matlab)
#' @author Carlos Camino (Ported version into R)
#'
#' @examples
#' spectral <- define.bands()
define.bands <- function() {
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



#' Compute the numerical Jacobian of the Fluspect leaf model with respect to its parameters
#'
#' @param x numeric vector. Leaf biochemistry/structure parameters at which to compute the numerical Jacobian.
#' @param spectral list. Spectral configuration used by the Fluspect model.
#' @param inputLeaf list. Leaf input parameters.
#' @param optipar list. Optical parameters used by the Fluspect model.
#'
#' @return J: a numeric array of partial derivatives (reflectance/transmittance with respect to each parameter in x).
#' @export
#'
#' @author 	Wout Verhoef, Christiaan van der Tol, Joris Timmermans (Original version in Matlab)
#' @author Carlos Camino (Ported version into R)
#'
#' @examples
#' \dontrun{
#' get.numjacobian(x, spectral, inputLeaf, optipar)
#' }
get.numjacobian <- function(x, spectral, inputLeaf, optipar) {
  stop(
    "get.numjacobian() is currently broken in three independent ways: ",
    "(1) it references 'leafbio', which is not one of its parameters (only 'inputLeaf' is) - ",
    "would fail with 'object leafbio not found'; ",
    "(2) it calls ToolsRTM::calc_fluspect_bcar(), which is itself broken (calls a function ",
    "that doesn't exist - see ToolsRTM's fluspect_functions.R); ",
    "(3) it calls ToolsRTM::get.calc_fluspect_bcar(), but that function is defined in THIS ",
    "package (SCOPEinR), not in ToolsRTM - wrong namespace. Needs review before use."
  )
  n <- length(x)
  res <- ToolsRTM::calc_fluspect_bcar(x, spectral, leafbio, optipar)
  fx <- c(res$refl, res$tran)
  step <- 1e-6
  J <- array(0, dim = c(length(fx), 2, n))

  for (i in 1:n) {
    xstep <- x
    xstep[i] <- x[i] + step
    res_step <- get.calc_fluspect_bcar(xstep, spectral, inputLeaf, optipar)
    fxstep <- c(res_step$refl, res_step$tran)
    J[, , i] <- (fxstep - fx) / step
  }

  return(J)
}

#' Run the combined B/Cx Fluspect leaf model from a packed parameter vector
#'
#' @param params numeric vector of length 6, in order: Cab, Cdm, Cw, Cs, Cca, N - leaf biochemistry/structure parameters.
#' @param spectral list. Spectral configuration used by the Fluspect model.
#' @param leafbio list. Baseline leaf parameters; Cab/Cdm/Cw/Cs/Cca/N are overwritten from params before simulation.
#' @param optipar list. Optical parameters used by the Fluspect model.
#'
#' @return A list with refl and tran (leaf reflectance/transmittance).
#' @export
#'
#' @examples
#' \dontrun{
#' get.calc_fluspect_bcar(params, spectral, leafbio, optipar)
#' }
get.calc_fluspect_bcar <- function(params, spectral, leafbio, optipar) {
  # Rest of the body removed as unreachable dead code (called a function,
  # ToolsRTM::fluspect_B_CX_PSI_PSII_combined(), that doesn't exist -- R CMD
  # check's static analysis flags the reference even though it's unreachable).
  stop(
    "get.calc_fluspect_bcar() is currently broken: it called ",
    "ToolsRTM::fluspect_B_CX_PSI_PSII_combined(), which does not exist in ToolsRTM. ",
    "The existing getFluspect.Cx()/getFluspect.B() (in ToolsRTM) use a different parameter ",
    "set (EWT, LMA, Cx...) than this function's leafbio (Cdm, Cw, Cca...), so they can't be ",
    "substituted without confirming the correct physical mapping between the two ",
    "parameterizations. Needs review before use - see project issue tracker / ask the package author."
  )
}




#' Cost function for fitting Fluspect leaf parameters to measured reflectance/transmittance
#'
#' \code{get.COST_4Fluspect.for.SCOPE} builds the residual vector between
#' simulated and measured leaf reflectance/transmittance for a given
#' parameter vector, intended for use inside a non-linear least-squares
#' optimizer when calibrating leaf biochemistry/structure parameters against
#' measured spectra. Note: this function calls an unqualified
#' \code{fluspect_B_CX_PSI_PSII_combined()}, which is neither defined in this
#' package nor exported by ToolsRTM - it will fail with "could not find
#' function" until that dependency is resolved. Not exported; needs review
#' before use.
#'
#' @param params numeric vector of length 8, in order: Cab, Cdm, Cw, Cs, Cca, Cant, Cx, N - candidate leaf biochemistry/structure parameters.
#' @param measurement list. Measured spectra, with elements \code{refl}, \code{tran} (measured reflectance/transmittance) and \code{std} (measurement uncertainty), all as a function of wavelength.
#' @param input list of length 6: \code{leafbio} (baseline leaf parameters), \code{optipar} (Fluspect optical parameters), \code{spectral} (spectral configuration, including \code{wlP} model wavelengths and \code{wlM} measurement wavelengths), \code{include} (list of 0/1 flags selecting which parameters in \code{params} are actually varied), \code{target} (character, \code{"0"} for combined reflectance+transmittance residuals, \code{"1"} for reflectance only, otherwise transmittance only), and \code{range} (list with \code{wlmin}/\code{wlmax} defining the wavelength range used to compute residuals).
#'
#' @return A list with: \code{er} (residual vector, or two-column matrix of reflectance/transmittance residuals when \code{target == "0"}), \code{refl}/\code{tran} (simulated reflectance/transmittance interpolated to \code{spectral$wlM}), and \code{leafopt} (full Fluspect output).
#'
#' @author Carlos Camino
#'
#' @examples
#' \dontrun{
#' cost <- get.COST_4Fluspect.for.SCOPE(params, measurement, input)
#' }
get.COST_4Fluspect.for.SCOPE <- function(params, measurement, input) {
  # Body removed as dead code -- called an unqualified fluspect_B_CX_PSI_PSII_combined(),
  # which is defined neither here nor in ToolsRTM (same root cause as
  # get.calc_fluspect_bcar() above). Previously had no stop() guard at all,
  # so calling this would fail with a raw "could not find function" instead
  # of an explanatory error -- fixed to fail clearly.
  stop(
    "get.COST_4Fluspect.for.SCOPE() is currently broken: it called an unqualified ",
    "fluspect_B_CX_PSI_PSII_combined(), which is defined neither in this package nor exported ",
    "by ToolsRTM. Needs review before use - see project issue tracker / ask the package author."
  )
}

