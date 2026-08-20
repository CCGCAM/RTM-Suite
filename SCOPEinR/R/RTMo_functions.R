
#' get.volscatt.scope version 2.0 from SCOPE model
#'
#' \code{volscatt.scope} calculates the scattering phase functions using a
#' radiative transfer model based on the optical properties of leaves. Specifically,
#' the model computes the fraction of radiation scattered in the forward and backward directions
#' for a given set of input parameters
#'
#' @param tts  Sun: zenith angle in degrees
#' @param tto  observation:zenith angle in degrees
#' @param psi   Difference of  azimuth angle between solar and viewing position
#' @param ttli leaf inclination array
#'
#' @return The function returns a list of four items: "chi_s", "chi_o", "frho", and "ftau".
#' These items represent the scattering phase functions for direct and diffuse radiation, respectively.
#'
#' @export
#' @details Date: 11 February 2008.
#' @author 	Wout Verhoef, Joris Timmermans (Original version in Matlab)
#' @author Carlos Camino (Ported version into R)
#'
#' @examples
#' get.volscatt.scope(tts = 30, tto = 10, psi = 0, ttli = seq(5, 85, by = 10))
get.volscatt.scope <- function(tts, tto, psi, ttli) {

  deg2rad <- pi/180
  nli     <- length(ttli)

  psi_rad         <- psi * deg2rad * rep(1, nli)

  cos_psi         <- cos(psi * deg2rad)                 # cosine of relative azimuth angle

  cos_ttli        <- cos(ttli * deg2rad)                # cosine of normal of upperside of leaf
  sin_ttli        <- sin(ttli * deg2rad)                # sine   of normal of upperside of leaf

  cos_tts         <- cos(tts * deg2rad)                 # cosine of sun zenith angle
  sin_tts         <- sin(tts * deg2rad)                 # sine   of sun zenith angle

  cos_tto         <- cos(tto * deg2rad)                 # cosine of observer zenith angle
  sin_tto         <- sin(tto * deg2rad)                 # sine   of observer zenith angle

  Cs              <- cos_ttli * cos_tts                 # p305{1}
  Ss              <- sin_ttli * sin_tts                 # p305{1}

  Co              <- cos_ttli * cos_tto                 # p305{1}
  So              <- sin_ttli * sin_tto                 # p305{1}

  As              <- pmax(Ss, Cs)
  Ao              <- pmax(So, Co)

  bts             <- acos(-Cs/As)                       # p305{1}
  bto             <- acos(-Co/Ao)                       # p305{2}

  chi_o           <- 2/pi * ((bto-pi/2) * Co + sin(bto) * So)
  chi_s           <- 2/pi * ((bts-pi/2) * Cs + sin(bts) * Ss)

  delta1          <- abs(bts-bto)                       # p308{1}
  delta2          <- pi - abs(bts + bto - pi)           # p308{1}

  Tot             <- psi_rad + delta1 + delta2           # pag 130{1}

  bt1             <- pmin(psi_rad, delta1)
  bt3             <- pmax(psi_rad, delta2)
  bt2             <- Tot - bt1 - bt3

  T1              <- 2*Cs*Co + Ss*So*cos_psi
  T2              <- sin(bt2) * (2*As*Ao + Ss*So*cos(bt1)*cos(bt3))

  Jmin            <- (bt2) * T1 - T2
  Jplus           <- (pi - bt2) * T1 + T2

  frho            <-  Jplus / (2 * pi^2)
  ftau            <- -Jmin / (2 * pi^2)

  # pag.309 wl-> pag 135{1}
  frho            <- pmax(0, frho)
  ftau            <- pmax(0, ftau)

  return(list(chi_s = chi_s, chi_o = chi_o, frho = frho, ftau = ftau))
}



#' get.Pso function
#'
#' \code{get.Pso} calculates the bi-directional gap probability \code{Pso}, i.e.
#' the joint probability of a sunlit and simultaneously viewed leaf (or soil
#' background) at a given normalized canopy depth, accounting for the
#' hot-spot effect between the solar and viewing directions.
#'
#' @param K numeric value. Canopy extinction coefficient in the direction of the observer.
#' @param k numeric value. Canopy extinction coefficient in the direction of the sun.
#' @param LAI numeric value. Total (one-sided) leaf area index of the canopy (m2 m-2).
#' @param q numeric value. Hot-spot size parameter (leaf width to canopy height ratio, dimensionless).
#' @param dso numeric value. Normalized distance between the sun and observer beams at the canopy level, derived from the solar/viewing zenith and relative azimuth angles (dimensionless).
#' @param xl numeric value. Normalized cumulative canopy depth at which \code{Pso} is evaluated (dimensionless, 0 at the top of the canopy to -1 at the bottom, expressed as a fraction of \code{LAI}).
#'
#' @return numeric value. Bi-directional gap probability \code{Pso} at depth \code{xl} (dimensionless, 0-1).
#' @export
#'
#' @author 	 Wout Verhoef, Christiaan van der Tol, Joris Timmermans (Original version in Matlab)
#' @author Carlos Camino (Ported version into R)
#'
#' @examples
#' \dontrun{
#' pso <- get.Pso(K, k, LAI, q, dso, xl)
#' }
get.Pso <- function(K, k, LAI, q, dso, xl) {
  if (dso != 0) {
    alf <- (dso/q) * 2/(k + K)
    pso <- exp((K + k) * LAI * xl + sqrt(K * k) * LAI/alf * (1 - exp(xl * alf)))
  } else {
    pso <- exp((K + k) * LAI * xl - sqrt(K * k) * LAI * xl)
  }
  return(pso)
}

#' get.reflectances
#'
#' \code{get.reflectances} propagates the thin-layer reflectance and
#' transmittance factors of the canopy (computed by the SAIL 4-stream
#' radiative transfer scheme) down through all \code{nl} canopy layers to the
#' soil, and back up, to obtain the directional-hemispherical and
#' hemispherical-hemispherical reflectance at the top of each layer.
#'
#' @param tau_ss numeric matrix with \code{nl} rows and \code{nwl} columns. Direct-direct transmittance of each thin canopy layer.
#' @param tau_sd numeric matrix with \code{nl} rows and \code{nwl} columns. Direct-diffuse transmittance of each thin canopy layer.
#' @param tau_dd numeric matrix with \code{nl} rows and \code{nwl} columns. Diffuse-diffuse transmittance of each thin canopy layer.
#' @param rho_dd numeric matrix with \code{nl} rows and \code{nwl} columns. Diffuse-diffuse reflectance of each thin canopy layer.
#' @param rho_sd numeric matrix with \code{nl} rows and \code{nwl} columns. Direct-diffuse reflectance of each thin canopy layer.
#' @param rsoil numeric vector of length \code{nwl}. Soil reflectance spectrum at the bottom boundary.
#' @param nl integer. Number of canopy layers.
#' @param nwl integer. Number of wavelengths in the spectral domain (default 2162, i.e. 400-2500 nm plus the thermal region).
#'
#' @return A list with:
#' \describe{
#'   \item{R_sd}{numeric matrix with \code{nl+1} rows and \code{nwl} columns. Directional-hemispherical reflectance at the top of each layer (including the soil, layer nl+1).}
#'   \item{R_dd}{numeric matrix with \code{nl+1} rows and \code{nwl} columns. Hemispherical-hemispherical reflectance at the top of each layer (including the soil, layer nl+1).}
#'   \item{Xss}{numeric vector of length \code{nl}. Direct-direct transmittance of each layer (identical to \code{tau_ss}, kept for use in flux propagation).}
#'   \item{Xsd}{numeric matrix with \code{nl} rows and \code{nwl} columns. Normalized direct-to-diffuse transmittance of each layer, accounting for multiple reflections with the layers below.}
#'   \item{Xdd}{numeric matrix with \code{nl} rows and \code{nwl} columns. Normalized diffuse-to-diffuse transmittance of each layer, accounting for multiple reflections with the layers below.}
#' }
#' @export
#'
#' @author 	 Wout Verhoef, Christiaan van der Tol, Joris Timmermans (Original version in Matlab)
#' @author Carlos Camino (Ported version into R)
#'
#' @examples
#' \dontrun{
#' refl <- get.reflectances(tau_ss, tau_sd, tau_dd, rho_dd, rho_sd, rsoil, nl, nwl = 2162)
#' }
get.reflectances <- function(tau_ss, tau_sd, tau_dd, rho_dd, rho_sd, rsoil, nl, nwl=2162) {


  if (nwl != 2162){
    ### check the wavelength for thermal region (>2500)
    message('please verify the length of the wavelengths (400:2500)')
    stop()


  }

    R_sd <- matrix(0, nl + 1, nwl)
    R_dd <- matrix(0, nl + 1, nwl)
    Xsd <- matrix(0, nl, nwl)
    Xdd <- matrix(0, nl, nwl)
    Xss <- rep(0, nl)
    # Initialize Xss as a matrix of zeros
    #Xss <- matrix(0, nrow = nl, ncol = length(k))

    R_sd[nl + 1, ] <- rsoil
    R_dd[nl + 1, ] <- rsoil


  #### Get RFL and Trans
  for (j in nl:1) {
    Xss[j] <- tau_ss[j]
    dnorm <- 1 - rho_dd[j, ] * R_dd[j + 1, ]
    Xsd[j, ] <- (tau_sd[j, ] + tau_ss[j] * R_sd[j + 1, ] * rho_dd[j, ]) / dnorm
    Xdd[j, ] <- tau_dd[j, ] / dnorm
    R_sd[j, ] <- rho_sd[j, ] + tau_dd[j, ] * (R_sd[j + 1, ] * Xss[j] + R_dd[j + 1, ] * Xsd[j, ])
    R_dd[j, ] <- rho_dd[j, ] + tau_dd[j, ] * R_dd[j + 1, ] * Xdd[j, ]
  }



  return(list(R_sd = R_sd, R_dd = R_dd, Xss = Xss, Xsd = Xsd, Xdd = Xdd))
}

#' get.fluxprofile
#'
#' \code{get.fluxprofile} propagates the top-of-canopy direct and diffuse
#' irradiance down through the \code{nl} canopy layers (and reflects them back
#' up), using the layer transmittance/reflectance factors from
#' \code{\link{get.reflectances}}, to obtain the vertical profile of direct
#' and diffuse radiation fluxes within the canopy.
#'
#' @param Esun_ numeric vector of length \code{nwl}. Top-of-canopy direct solar irradiance spectrum.
#' @param Esky_ numeric vector of length \code{nwl}. Top-of-canopy diffuse sky irradiance spectrum.
#' @param rsoil numeric vector of length \code{nwl}. Soil reflectance spectrum (bottom boundary condition).
#' @param Xss numeric vector of length \code{nl}. Direct-direct transmittance of each layer, as returned by \code{\link{get.reflectances}}.
#' @param Xsd numeric matrix with \code{nl} rows and \code{nwl} columns. Normalized direct-to-diffuse transmittance of each layer, as returned by \code{\link{get.reflectances}}.
#' @param Xdd numeric matrix with \code{nl} rows and \code{nwl} columns. Normalized diffuse-to-diffuse transmittance of each layer, as returned by \code{\link{get.reflectances}}.
#' @param R_sd numeric matrix with \code{nl+1} rows and \code{nwl} columns. Directional-hemispherical reflectance at the top of each layer, as returned by \code{\link{get.reflectances}}.
#' @param R_dd numeric matrix with \code{nl+1} rows and \code{nwl} columns. Hemispherical-hemispherical reflectance at the top of each layer, as returned by \code{\link{get.reflectances}}.
#' @param nl integer. Number of canopy layers.
#' @param nwl integer. Number of wavelengths in the spectral domain.
#' @param rs.thermal numeric. Soil reflectance used to initialize the flux matrices and (when \code{Xsd} does not cover the thermal region) to extend the spectral domain into the thermal region; default 0.06.
#'
#' @return A list with:
#' \describe{
#'   \item{Es_}{numeric matrix with \code{nl+1} rows and \code{nwl} columns. Direct solar irradiance at the top of each layer (down to the soil).}
#'   \item{Emin_}{numeric matrix with \code{nl+1} rows and \code{nwl} columns. Downward diffuse irradiance at the top of each layer.}
#'   \item{Eplu_}{numeric matrix with \code{nl+1} rows and \code{nwl} columns. Upward diffuse irradiance at the top of each layer.}
#' }
#' @export
#'
#' @author 	 Wout Verhoef, Christiaan van der Tol, Joris Timmermans (Original version in Matlab)
#' @author Carlos Camino (Ported version into R)
#'
#' @examples
#' \dontrun{
#' Eflux <- get.fluxprofile(Esun_, Esky_, rsoil, Xss, Xsd, Xdd, R_sd, R_dd, nl, nwl)
#' }
get.fluxprofile <- function(Esun_, Esky_, rsoil, Xss, Xsd, Xdd, R_sd, R_dd, nl, nwl, rs.thermal=0.06) {

  if (missing(rs.thermal)){
    rs_thermal = 0.01
  } else{
    rs_thermal = rs.thermal
  }

  #nl <- nl + 1 # increment nl to match indexing
  Es_ <- matrix(rs_thermal, nrow = nl+1, ncol = nwl)
  Emin_ <- matrix(rs_thermal, nrow = nl+1, ncol = nwl)
  Eplu_ <- matrix(rs_thermal, nrow = nl+1, ncol = nwl)

  ### adding value of
  if (dim(Xsd)[2] == 2001){


    rs_vec <- 2162 - 2001
    r_thermal <-matrix(rs_thermal, ncol = rs_vec, nrow= nrow(R_dd))
    r_thermal.X <-matrix(rs_thermal, ncol = rs_vec, nrow= nrow(Xsd))

    ### convert rfl to a matrix with columns ==nwl (all spectrla domain addinfg thermal)
    # repeat rs_thermal to match the number of columns to add
    Xsd <- cbind(Xsd, r_thermal.X)
    Xdd <- cbind(Xdd, r_thermal.X)
    R_dd <- cbind(R_dd, r_thermal)
    R_sd <- cbind(R_sd, r_thermal)
    rsoil = c(rsoil, rep(rs_thermal,rs_vec))
  }


  Es_[1,] <- Esun_
  Emin_[1,] <- Esky_

  for (j in c(1:nl)) {

    Es_[j+1,] <- Xss * Es_[j,] #matrix(Xss,ncol=1)[j] %*% (Es_[j,]) --> before
    #Emin_[j+1,] <- Xsd[j,] * Es_[j,] + matrix(Xdd,ncol=1)[j] %*% Emin_[j,] --> before
    Emin_[j+1,] <- Xsd[j,] * Es_[j,] + Xdd[j,] * Emin_[j,]
    Eplu_[j,] <- R_sd[j,] * Es_[j,] + R_dd[j,] * Emin_[j,]

  }
  #Eplu_[nl,] <- rsoil * (Es_[nl,] + Emin_[nl,])
  Eplu_[nl + 1, ] <- rsoil * (Es_[nl + 1, ] + Emin_[nl + 1, ])
  # CvdT added calculation of Eplu_[nl,]

  return(list(Es_ = Es_, Emin_ = Emin_, Eplu_ = Eplu_))
}


#' get.calcTOCirr
#'
#' \code{get.calcTOCirr} calculates the top-of-canopy (TOC) direct solar
#' (\code{Esun_}) and diffuse sky (\code{Esky_}) irradiance spectra. If
#' \code{atmo} already contains pre-computed \code{Esun_}/\code{Esky_} these
#' are returned directly; otherwise they are derived from the MODTRAN
#' atmospheric transmittance/radiance functions (T1-T16) evaluated at the
#' SCOPE wavelengths, optionally rescaled to match the observed broadband
#' incoming shortwave (\code{Rin}) and longwave (\code{Rli}) radiation.
#'
#' @param atmo list or data.frame. Atmospheric input, either pre-computed \code{Esun_}/\code{Esky_} spectra, or MODTRAN output containing wavenumber (\code{WN}) and transmittance/radiance functions \code{T1, T2, T3, T4, T5, T12, T16}.
#' @param meteo list. Meteorological data; must contain \code{Ta} (air temperature, deg C) and, when rescaling is required, \code{Rin} (incoming shortwave radiation, W m-2) and \code{Rli} (incoming longwave radiation, W m-2). Use \code{-999} for \code{Rin} to skip rescaling.
#' @param rdd numeric vector of length \code{nwl}. Top-of-canopy hemispherical-hemispherical reflectance spectrum.
#' @param rsd numeric vector of length \code{nwl}. Top-of-canopy directional-hemispherical reflectance spectrum.
#' @param wl numeric vector of length \code{nwl}. Wavelengths of the spectral domain (nm).
#' @param nwl integer. Number of wavelengths in the spectral domain.
#'
#' @return A list with:
#' \describe{
#'   \item{Esun_}{numeric vector of length \code{nwl}. Top-of-canopy direct solar irradiance spectrum.}
#'   \item{Esky_}{numeric vector of length \code{nwl}. Top-of-canopy diffuse sky irradiance spectrum.}
#' }
#' @export
#'
#'
#' @author 	 Wout Verhoef, Christiaan van der Tol, Joris Timmermans (Original version in Matlab)
#' @author Carlos Camino (Ported version into R)
#'
#' @examples
#' \dontrun{
#' TOCirr <- get.calcTOCirr(atmo, meteo, rdd, rsd, wl, nwl)
#' }
get.calcTOCirr <- function(atmo, meteo, rdd, rsd, wl, nwl){

  Fd      <- rep(0, nwl)
  Ls      <- get.Planck(wl, meteo[['Ta']] + 273.15)

  ##### This function need to be revised with a MODtrans file


  if(!"Esun_" %in% names(atmo)){

    wln_ <-c(spectral$reg1,spectral$reg2,spectral$reg3)
    atmo.selected <- atmo[atmo$WN %in% wl, c('WN','T1', 'T2','T3','T4','T5','T12','T16')]

    t1  <- atmo.selected[, 'T1']
    t3  <- atmo.selected[, 'T3']
    t4  <- atmo.selected[, 'T4']
    t5  <- atmo.selected[, 'T5']
    t12 <- atmo.selected[, 'T12']
    t16 <- atmo.selected[, 'T16']

    # radiation fluxes, downward and upward (these all have dimenstion [nwl]
    # first calculate hemispherical reflectances rsd and rdd according to SAIL
    # these are assumed for the reflectance of the surroundings
    # rdo is computed with SAIL as well
    # assume Fd of surroundings = 0 for the momemnt
    # initial guess of temperature of surroundings from Ta;

    Esun_   <- pmax(1e-6, pi*t1*t4)

    rdd_wl <- c(c(rdd), rep(soil[['rs_thermal']],length(spectral[['IwlT']])))
    rsd_wl <- c(c(rsd), rep(soil[['rs_thermal']],length(spectral[['IwlT']])))
    Esky_   <- pmax(1e-6, pi/(1-t3*rdd_wl)*(t1*(t5+t12*rsd_wl)+Fd+(1-rdd_wl)*Ls*t3+t16))

    # fractional contributions of Esun and Esky to total incident radiation in
    # optical and thermal parts of the spectrum
    if(meteo[['Rin']] != -999){
      # fractional contributions of Esun and Esky to total incident radiation in
      # optical and thermal parts of the spectrum

      fEsuno  <- rep(0, nwl)
      fEskyo  <- rep(0, nwl)
      fEsunt  <- rep(0, nwl)
      fEskyt  <- rep(0, nwl)

      J_o             <- which(wl<3000)             #find optical spectrum
      Esunto          <- 0.001 * Sint(Esun_[J_o], wl[J_o])  #Calculate optical sun fluxes (by Integration), including conversion mW >> W
      Eskyto          <- 0.001 * Sint(Esky_[J_o], wl[J_o])  #Calculate optical sun fluxes (by Integration)
      Etoto           <- Esunto + Eskyto             #Calculate total fluxes
      fEsuno[J_o]     <- Esun_[J_o]/Etoto            #fraction of contribution of Sun fluxes to total light
      fEskyo[J_o]     <- Esky_[J_o]/Etoto            #fraction of contribution of Sky fluxes to total light

      J_t             <- which(wl>=3000)            #find thermal spectrum
      Esuntt          <- 0.001 * Sint(Esun_[J_t], wl[J_t])  #Themal solar fluxes
      Eskytt          <- 0.001 * Sint(Esky_[J_t], wl[J_t])  #Thermal Sky fluxes
      Etott           <- Eskytt + Esuntt             #Total
      fEsunt[J_t]     <- Esun_[J_t]/Etott            #fraction from Esun
      fEskyt[J_t]     <- Esky[J_t]/Etott  #fraction from Esky


      Esun_[J_o] = fEsuno[J_o] * meteo[['Rin']]
      Esky_[J_o] = fEskyo[J_o] * meteo[['Rin']]
      Esun_[J_t] = fEsunt[J_t] * meteo[['Rli']]
      Esky_[J_t] = fEskyt[J_t] * meteo[['Rli']]

    }
  } else {
    Esun_ = atmo[['Esun_']]
    Esky_ = atmo[['Esky_']]
  }

  return(list(Esun_ = Esun_, Esky_ = Esky_))
}


