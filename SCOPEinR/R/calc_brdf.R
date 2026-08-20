
#' @title get.brdf
#'
#' @description
#' \code{get.brdf} simulates the bidirectional reflectance (and, optionally,
#' the directional thermal, fluorescence and xanthophyll-related radiance) of
#' the canopy over a large set of viewing angles, by repeatedly calling
#' \code{\link{getRTMo}} (and, depending on \code{data.opts}, \code{get.RTMt.planck},
#' \code{get.RTMf}, \code{get.RTMz}) for each requested angle combination. The
#' viewing angles requested in \code{data.directional} are extended with a set
#' of angles for hot-spot oversampling and principal-plane oversampling, and
#' duplicate angle combinations are removed before the loop.
#'
#' @param data.spectral list. Spectral band definitions used throughout the model (e.g. \code{wlS}, \code{wlF}, \code{wlT}).
#' @param data.angles list. Base observation geometry, with element \code{tts} (solar zenith angle, degrees) used to construct the hot-spot oversampling angles.
#' @param data.rad list. Radiation fluxes and top-of-canopy irradiance computed previously by \code{getRTMo}, passed through to the thermal/fluorescence/xanthophyll sub-models.
#' @param data.directional list. User-requested directional observation geometry, with elements \code{psi} (relative azimuth angles, degrees) and \code{tto} (viewing zenith angles, degrees).
#' @param atmo list or data.frame. Atmospheric input (MODTRAN transmittance functions or precomputed \code{Esun_}/\code{Esky_}), passed to \code{getRTMo}.
#' @param data.soil list. Soil reflectance and related soil properties.
#' @param data.leafopt list. Leaf optical properties (reflectance/transmittance spectra) computed by the leaf optical model.
#' @param data.leafbio list. Leaf biochemical/biophysical parameters (\code{Cab}, \code{Cw}, \code{Cdm}, \code{N}, ...).
#' @param data.canopy list. Canopy structural properties (LAI, leaf inclination distribution, hot-spot parameter, etc.).
#' @param data.gap list. Gap fraction probabilities (\code{Ps}, \code{Po}, \code{Pso}) needed to combine sunlit/shaded contributions in the thermal, fluorescence and xanthophyll sub-models.
#' @param data.meteo list. Meteorological forcing data (e.g. air temperature, vapour pressure).
#' @param data.thermal list. Component temperatures from the energy balance: \code{Tcu} (sunlit canopy), \code{Tch} (shaded canopy), \code{Tsu} (sunlit soil), \code{Tsh} (shaded soil), all in deg C.
#' @param data.bcu list. Biochemical outputs for sunlit leaves, including \code{eta} (relative fluorescence emission efficiency) and \code{Kn} (non-photochemical quenching rate constant).
#' @param data.bch list. Biochemical outputs for shaded leaves, with the same structure as \code{data.bcu}.
#' @param data.opts data.frame. Model configuration options; rows used here are \code{calc_planck} (compute thermal radiance spectrum), \code{calc_fluor} (compute chlorophyll fluorescence), and \code{calc_xanthophyllabs} (include xanthophyll de-epoxidation effect on reflectance).
#' @param get.plots logical. If \code{TRUE}, diagnostic plots are produced by the underlying RTM calls. Default \code{FALSE}.
#'
#' @return A list \code{directional} with matrices indexed by wavelength (rows) and viewing angle (columns), including:
#' \describe{
#'   \item{psi, tto}{numeric vectors. The (deduplicated) relative azimuth and viewing zenith angles actually simulated.}
#'   \item{refl_}{TOC reflectance spectrum per angle.}
#'   \item{rso_}{Directional-directional BRDF per angle.}
#'   \item{Lo_}{TOC outgoing radiance per angle.}
#'   \item{LoF_}{TOC fluorescence radiance per angle (if \code{calc_fluor} is enabled).}
#'   \item{Lot_}{TOC thermal radiance per angle (if \code{calc_planck} is enabled).}
#'   \item{Eoutte, BrightnessT}{Placeholder outputs (currently left at their initialized zero values).}
#' }
#' @export
#'
#' @author Carlos Camino
#'
#' @examples
#' \dontrun{
#' directional <- get.brdf(data.spectral, data.angles, data.rad, data.directional, atmo,
#'                          data.soil, data.leafopt, data.leafbio, data.canopy, data.gap,
#'                          data.meteo, data.thermal, data.bcu, data.bch, data.opts)
#' }
get.brdf <- function(data.spectral,data.angles,data.rad,
                     data.directional,
                     atmo,
                     data.soil,data.leafopt,data.leafbio,data.canopy,data.gap,
                     data.meteo,data.thermal,
                     data.bcu,data.bch,data.opts,
                     get.plots=F) {




  # simulates observations from a large number of viewing angles
  # modified: 30 April 2020, CvdT, removed repeated angle combinations.

  if (missing(data.opts)){
    stop('please use options for Lite, Calc_vert_profiles ')

  } else{

    options.calc_planck         = data.opts[3,]   # calculate spectrum of thermal radiation
    options.calc_fluor          = data.opts[2,]   # calculate chlorophyll fluorescence in observation direction
    options.calc_xanthophyllabs = data.opts[4,]   # include simulation of reflectance dependence on de-epoxydation state

  }



  if (missing(get.plots)){
    get.plots = FALSE
  }
  ## input
  tts <- data.angles$tts
  psi_hot <- c(0, 0, 0, 0, 0, 2, 358) # [noa_o] angles for hotspot oversampling
  tto_hot <- c(tts, tts + 2, tts + 4, tts - 2, tts - 4, tts, tts) # [noa_o] angles for hotspot oversampling

  psi_plane <- c(rep(0, 6), rep(180, 6), rep(90, 6), rep(270, 6)) # angles for plane oversampling
  tto_plane <- c(10:60, 10:60, 10:60, 10:60) # angles for plane oversampling

  psi <- c(data.directional[['psi']], psi_hot, psi_plane)
  tto <- c(data.directional[['tto']] , tto_hot, tto_plane)
  directional <- list()
  ## remove duplicates
  unique_angles <- unique(cbind(psi, tto))
  directional[['psi']] <- unique_angles[,1]
  directional[['tto']] <- unique_angles[,2]
  na <- length(unique_angles[,1])

  ## allocate memory
  directional[['brdf_']] <- matrix(0, nrow=length(data.spectral$wlS), ncol=na) # [nwlS, no of angles]
  directional[['Eoutte']] <- matrix(0, nrow=1, ncol=na) # [1, no of angles]
  directional[['BrightnessT']] <- matrix(0, nrow=1, ncol=na) # [1, no of angles]

  directional[['LoF_']] <- matrix(0, nrow=length(data.spectral$wlF), ncol=na) # [nwlF, no of angles]
  directional[['Lot_']] <- matrix(0, nrow=length(data.spectral$wlT), ncol=na) # [nwlF, no of angles]

  directional[['refl_']] <- matrix(0, nrow=length(data.spectral$wlS), ncol=na) # [nwlF, no of angles]
  directional[['rso_']] <- matrix(0, nrow=length(data.spectral$wlS), ncol=na)
  directional[['Lo_']] <- matrix(0, nrow=length(data.spectral$wlS), ncol=na)
  ## other preparations
  directional_angles <- data.angles

  ## loop over the angles

  progress_bar = txtProgressBar(min=0, max=na, style = 3, char="=")

  for (j in 1:na) {

    setTxtProgressBar(progress_bar, j)
    # optical BRDF
    directional_angles[['tto']] <- directional[['tto']][j]
    directional_angles[['psi']] <- directional[['psi']][j]

    directional.getRTM0<-getRTMo(data.spectral,atmo,data.soil,data.leafopt,data.canopy,data.leafbio,
                       data.angles=directional_angles,data.meteo,data.opts=data.opts,get.plots=F)

    directional.rad <-directional.getRTM0$data.rad


    directional[['refl_']][,j] <- directional.rad[['refl']] # [nwl] reflectance (spectral) (nm-1)
    directional[['rso_']][,j] <- directional.rad[['rso']] # [nwl] BRDF (spectral) (nm-1)

    # thermal directional brightness temperatures (Planck)
    if (options.calc_planck$Value == 1) {


      directional.RTMt.planck <- get.RTMt.planck(data.spectral=data.spectral,data.rad=data.rad,data.soil=data.soil,
                                                 data.leafbio=data.leafbio, data.leafopt=data.leafopt,
                                                 data.canopy=data.canopy,
                                                 data.gap=data.gap,Tcu=data.thermal[['Tcu']],Tch=data.thermal[['Tch']],
                                  Tsu=data.thermal[['Tsu']],Tsh=data.thermal[['Tsh']],
                                  get.plots=F)

      directional[['Lot_']][,j] <- directional.RTMt.planck$Lot_[data.spectral$IwlT] + directional.RTMt.planck$Lo_[data.spectral$IwlT] # [nwlt] emitted plus reflected diffuse radiance at

    }

    if (options.calc_fluor$Value == 1) {

      directional.RTMf <- get.RTMf(data.spectral=data.spectral,data.rad=data.rad,
                                   data.soil=data.soil,data.leafopt=data.leafopt,
                                   data.canopy=data.canopy,data.gap=data.gap,data.angles=data.angles,
                           #relative fluorescence emission efficiency for sunlit leaves
                           data.etau=data.bcu[['eta']],
                           #relative fluorescence emission efficiency for shaded leaves
                           data.etah =data.bch[['eta']],get.plots=F)


      directional$LoF_[, j] <- directional.RTMf$LoF_
    }

    if (options.calc_xanthophyllabs$Value == 1) {

      directional.RTMz <- get.RTMz(data.spectral=data.spectral,data.rad=data.rad,
                                   data.soil=data.soil,data.leafopt=data.leafopt,
                                   data.canopy=data.canopy,data.gap=data.gap,data.angles=data.angles,
                           data.Knu = data.bcu[['Kn']],data.Knh = data.bch[['Kn']],
                           get.plots=F)



     directional$Lo_[,j] <- directional.RTMz$Lo_

    }
  } # end loop na for all angles
  close(progress_bar)
  return(directional)

} # end angles
