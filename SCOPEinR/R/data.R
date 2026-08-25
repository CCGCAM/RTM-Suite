#' Physical and physiological constants used by the SCOPE model
#'
#' A lookup table of physical constants (e.g. molar masses, gas constant,
#' Stefan-Boltzmann constant) used throughout the energy balance and
#' radiative transfer routines.
#'
#' @format A data frame with one row per constant, columns \code{constant} (name) and a value column.
#' @source Van der Tol, C., W. Verhoef, J. Timmermans, A. Verhoef, and Z. Su. 2009. Biogeosciences 6(12): 3109-29.
"constants"

#' Default SCOPE run options
#'
#' Default values for the SCOPE model run-time options (equivalent to \code{setoptions.csv}),
#' controlling which radiative transfer/energy balance/output modules are enabled.
#'
#' @format A data frame with one row per option and a \code{Value} column.
#' @source Yang, P., E. Prikaziuk, W. Verhoef, and C. van der Tol. 2020. Geoscientific Model Development Discussions.
"data.opts"

#' Leaf optical parameters (PROSPECT/Fluspect, various parameterizations)
#'
#' Wavelength-dependent leaf-level optical coefficients (refractive index,
#' specific absorption coefficients per pigment, PSI/PSII quantum yield
#' spectra) required by the Fluspect/PROSPECT leaf models. Several versions
#' are bundled, corresponding to different published parameterizations.
#'
#' @format A data frame/list with one row per wavelength and one column per optical coefficient.
#' @source \code{optipar}, \code{optipar.2015}: Feret et al. 2008 (PROSPECT-4); \code{optipar2017.ProspectD}: Feret et al. 2017 (PROSPECT-D);
#'   \code{optipar2020.prospectD.BSM2019}: Feret et al. 2017 + BSM soil model (Yang et al. 2019); \code{optipar2021.Pro.CX}: PROSPECT-PRO (Feret et al. 2021).
#' @name optipar-datasets
#' @aliases optipar optipar.2015 optipar2017.ProspectD optipar2020.prospectD.BSM2019 optipar2021.Pro.CX
"optipar"

#' @rdname optipar-datasets
"optipar.2015"

#' @rdname optipar-datasets
"optipar2017.ProspectD"

#' @rdname optipar-datasets
"optipar2020.prospectD.BSM2019"

#' @rdname optipar-datasets
"optipar2021.Pro.CX"

#' Example meteorological/radiation time series for diurnal SCOPE simulations
#'
#' Example time series (one value per time step) used by \code{\link{getLUT_time}}
#' to drive multi-timestep SCOPE runs: incoming shortwave radiation
#' (\code{Rin_}), incoming longwave radiation (\code{Rli_}), sun/sky
#' irradiance spectra (\code{Esun_}/\code{Esky_}), air temperature
#' (\code{Ta_}), vapour pressure (\code{ea_}), air pressure (\code{p_}),
#' wind speed (\code{u_}), time-of-day (\code{t_}) and year (\code{year_}).
#'
#' @format A numeric vector or data frame, one row/element per time step.
#' @name meteo-timeseries
#' @aliases Rin_ Rli_ Esun_ Esky_ Ta_ ea_ p_ u_ t_ year_
"Rin_"

#' @rdname meteo-timeseries
"Rli_"

#' @rdname meteo-timeseries
"Esun_"

#' @rdname meteo-timeseries
"Esky_"

#' @rdname meteo-timeseries
"Ta_"

#' @rdname meteo-timeseries
"ea_"

#' @rdname meteo-timeseries
"p_"

#' @rdname meteo-timeseries
"u_"

#' @rdname meteo-timeseries
"t_"

#' @rdname meteo-timeseries
"year_"

#' BRDF viewing/illumination angle configurations
#'
#' Predefined sets of solar/viewing zenith and relative azimuth angles used
#' to compute the bidirectional reflectance distribution function (BRDF) of
#' the canopy. \code{brdf_angles_no_oversampling} omits the finer angular
#' oversampling used by \code{brdf_angles}/\code{brdf_angles2}.
#'
#' @format A data frame with columns for solar zenith, viewing zenith and relative azimuth angles (degrees).
#' @name brdf-angle-datasets
#' @aliases brdf_angles brdf_angles2 brdf_angles_no_oversampling
"brdf_angles"

#' @rdname brdf-angle-datasets
"brdf_angles2"

#' @rdname brdf-angle-datasets
"brdf_angles_no_oversampling"

#' Default SCOPE input look-up table
#'
#' A ready-to-use look-up table of SCOPE model inputs (leaf biochemistry,
#' canopy structure, viewing/illumination geometry, meteorology) for a
#' single default simulation, used in examples and the package tutorials.
#'
#' @format A data frame, one row per simulation, one column per input parameter.
"SCOPE.LUT.default"

#' Example SCOPE input LUT and input border/range definitions
#'
#' \code{inputsSCOPE} is an example base look-up table of SCOPE input
#' parameters (as read from \code{inst/input/inputs_SCOPE.csv}, used by
#' \code{\link{getLUT.SCOPE}} to generate randomized LUTs). \code{input_border}
#' gives the valid minimum/maximum range for each input parameter.
#'
#' @format A data frame, one row per input parameter, with default value and (for \code{input_border}) min/max range columns.
#' @name scope-input-datasets
#' @aliases inputsSCOPE input_border
"inputsSCOPE"

#' @rdname scope-input-datasets
"input_border"

#' Example leaf reflectance/transmittance spectrum
#'
#' An example leaf optical spectrum (reflectance and/or transmittance vs.
#' wavelength) used for testing and illustrating the leaf optical models.
#'
#' @format A data frame with a wavelength column and one or more reflectance/transmittance columns.
"leaf_spectrum"

#' Example soil reflectance spectrum
#'
#' An example dry/wet soil reflectance spectrum (400-2400 nm), as used by
#' the BSM soil model and the SCOPE radiative transfer routines.
#'
#' @format A data frame or vector of reflectance values, one per wavelength.
"soil.rfl"

#' Biochemical model lookup tables for Vcmax25/Jmax temperature response
#'
#' Lookup tables of the temperature-response parameters used by the
#' biochemical (photosynthesis) model to scale maximum carboxylation
#' (\code{Vcmax25}) and electron transport (\code{Jmax}) capacity with leaf
#' temperature.
#'
#' @format A data frame of temperature-response coefficients.
#' @name biochemistry-tables
#' @aliases table_Vcmax_ table_Jmax_
"table_Vcmax_"

#' @rdname biochemistry-tables
"table_Jmax_"
