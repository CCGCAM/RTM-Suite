#' PROSPECT-family leaf optical specific absorption coefficients (PDB set)
#'
#' Wavelength-indexed specific absorption/refractive-index coefficients used
#' by the PROSPECT-PRO / PROSPECT-D leaf optical model implementations in
#' this package (see [prospect_DB()]).
#'
#' @format A data frame with 2101 rows (400-2500 nm, 1 nm steps) and 12
#'   columns: `wavelength`, `Refrac_leafm` (leaf refractive index), specific
#'   absorption coefficients for chlorophyll (`SC_chl`), carotenoids
#'   (`SC_car`), anthocyanins (`SC_Anth`), brown pigments (`SC_Brwon`),
#'   water (`SC_Cw`) and dry matter (`SC_Cm`), plus direct/diffuse solar
#'   irradiance and dry/wet soil reflectance reference spectra.
#' @keywords datasets
"dataSpec_PDB"

#' PROSPECT leaf optical specific absorption coefficients (PRO set)
#'
#' Alternative wavelength-indexed specific absorption coefficient set for
#' the PROSPECT-PRO leaf model, following Feret et al.'s protein/non-protein
#' dry matter decomposition.
#'
#' @format A data frame with 2101 rows (400-2500 nm, 1 nm steps) and 10
#'   columns: `wave`, leaf refractive index (`RI_leafmaterial`), and
#'   specific absorption coefficients for chlorophyll, carotenoids,
#'   anthocyanins, brown pigments, water, dry matter, protein and
#'   non-protein dry matter fractions.
#' @keywords datasets
"dataSpec_PRO"

#' Fluspect specific absorption coefficients (2003/2008 parameterisation)
#'
#' Leaf optical coefficient table used by the Fluspect fluorescence-emission
#' model implementation, including chlorophyll fluorescence excitation
#' efficiency spectra.
#'
#' @format A data frame with 2001 rows (400-2400 nm) and 20 columns
#'   including wavelength (`wl`), refractive index (`nr`), specific
#'   absorption coefficients (`Kab`, `Kca`, `Ks`, `Kw`, `Kdm`), fluorescence
#'   excitation efficiencies (`phiI`, `phiII`), carotenoid violaxanthin/
#'   zeaxanthin coefficients (`KcaV`, `KcaZ`), anthocyanin coefficient
#'   (`Kant`) and generalised soil/senescence vectors (`GSV.1`, `GSV.2`).
#' @keywords datasets
"optipar"

#' Fluspect specific absorption coefficients (2015 parameterisation)
#'
#' Updated Fluspect coefficient table (Vilfan et al. 2016 parameterisation)
#' used by the mSCOPE/Fluspect-CX fluorescence model path.
#'
#' @format A data frame with 2001 rows (400-2400 nm) and 15 columns
#'   including wavelength (`wl`), refractive index (`nr`), specific
#'   absorption coefficients (`Kab`, `Kcar`, `Ks`, `Kw`, `Kdm`, `nw`),
#'   fluorescence basis spectra (`flu1`, `flu2`), three generalised soil
#'   vectors (`GSV1`-`GSV3`) and violaxanthin/zeaxanthin carotenoid
#'   coefficients (`KcaV`, `KcaZ`).
#' @keywords datasets
"optipar.2015"

#' Liberty leaf model specific absorption coefficients
#'
#' Wavelength-indexed coefficient table used by the Liberty needle/leaf
#' optical model implementation (see `liberty()`).
#'
#' @format A data frame with 421 rows and 5 unnamed numeric columns
#'   (`V1`-`V5`) holding the wavelength grid and the model's specific
#'   absorption/scattering coefficients, as distributed with the original
#'   Liberty model source data.
#' @keywords datasets
"dataspec.liberty"

#' Example measured leaf reflectance/transmittance spectrum
#'
#' A single example leaf optical spectrum bundled for demonstrating spectral
#' plotting and convolution functions in the package vignettes.
#'
#' @format A data frame with 2100 rows (one per 1 nm wavelength band,
#'   400-2500 nm) and 3 columns holding wavelength and paired
#'   reflectance/transmittance-like values.
#' @keywords datasets
"leaf_spectrum"

#' RTM parameter ranges and prior distributions (leaf/canopy models, shared)
#'
#' `inputsFlUSPECT`, `inputsINFORM`, `inputsLiberty`, `inputsPROSAIL`,
#' `inputsSPART` and `inputsRTMs` each give the default parameter bounds and
#' sampling distributions used by [getLUT()] and related look-up-table
#' generators for the corresponding model (Fluspect, INFORM, Liberty,
#' PROSAIL, SPART, and a combined table across models respectively).
#'
#' @format A data frame, one row per model parameter, with columns
#'   `variable` (parameter name), `lower`/`upper` (sampling bounds),
#'   `units`, `Distribution` (sampling distribution family),
#'   `Mean_D`/`Std_D` (distribution mean/SD where applicable),
#'   `Dependencies` (other parameters this one covaries with, if any),
#'   `use.default` (whether the package default is used) and `default`
#'   (default value). `inputsRTMs` additionally has a `model` column
#'   identifying which RTM each row belongs to.
#' @keywords datasets
#' @name rtm-input-specs
"inputsFlUSPECT"

#' @rdname rtm-input-specs
"inputsINFORM"

#' @rdname rtm-input-specs
"inputsLiberty"

#' @rdname rtm-input-specs
"inputsPROSAIL"

#' @rdname rtm-input-specs
"inputsSPART"

#' @rdname rtm-input-specs
"inputsRTMs"

#' Satellite sensor band, SMAC atmospheric-correction and SRF metadata
#'
#' `LANDSAT4.TM`, `LANDSAT5.TM`, `LANDSAT7.ETM`, `LANDSAT8.OLI`,
#' `Sentinel2A.MSI`, `Sentinel2B.MSI`, `Sentinel3A.OLCI`, `Sentinel3B.OLCI`
#' and `TerraAqua.MODIS` each bundle the per-band metadata, SMAC atmospheric
#' correction coefficients and spectral response function (SRF) tables
#' needed by [get.spectral.convolution()] and the sensor-convolution
#' vignettes to simulate at-sensor reflectance for that mission.
#'
#' @format A named list with entries that include `mission`/`name`
#'   (sensor identifiers), `band_id_all` (band names), `res_spatials`
#'   (spatial resolution per band), `rang_wvls` (nominal spectral range per
#'   band), `swath_widths`, `revisit_days`/`revisit_time`, `band_width`,
#'   `center_wvl` (band centre wavelengths), `SMAC_coef` (SMAC
#'   atmospheric-correction coefficient table), `wl_smac`/`p_srf`/`wl_srf`
#'   (or the `wvl_srf`/`p_srf_smac`/`wl_srf_smac` equivalents for the
#'   Sentinel-2 entries) giving the per-band spectral response function
#'   sampled on the SMAC/native wavelength grids, and `id_smac_in_all`
#'   mapping SMAC bands to the full band list. Exact component names vary
#'   slightly by sensor family; see the source list names for the specific
#'   object.
#' @keywords datasets
#' @name sensor-smac-srf
"LANDSAT4.TM"

#' @rdname sensor-smac-srf
"LANDSAT5.TM"

#' @rdname sensor-smac-srf
"LANDSAT7.ETM"

#' @rdname sensor-smac-srf
"LANDSAT8.OLI"

#' @rdname sensor-smac-srf
"Sentinel2A.MSI"

#' @rdname sensor-smac-srf
"Sentinel2B.MSI"

#' @rdname sensor-smac-srf
"Sentinel3A.OLCI"

#' @rdname sensor-smac-srf
"Sentinel3B.OLCI"

#' @rdname sensor-smac-srf
"TerraAqua.MODIS"

#' EnMAP hyperspectral sensor band characteristics
#'
#' Per-channel centre wavelength and full-width-half-maximum (FWHM) table
#' for the EnMAP hyperspectral imager, used by the sensor-convolution
#' functions to build a Gaussian spectral response for EnMAP bands.
#'
#' @format A data frame with 242 rows (one per EnMAP channel) and 4 columns:
#'   `Sensor`, `channel` (band index), `center` (centre wavelength, nm) and
#'   `fwhm` (full width at half maximum, nm).
#' @keywords datasets
"EnMap.characteristics"

#' PRISMA hyperspectral sensor band FWHM table
#'
#' Per-channel centre wavelength and full-width-half-maximum (FWHM) table
#' for the PRISMA hyperspectral imager, used to build a Gaussian spectral
#' response for PRISMA bands.
#'
#' @format A data frame with 234 rows (one per PRISMA channel) and 3
#'   columns: `QB` (band quality/type flag), `wavelength` (centre
#'   wavelength, nm) and `fwhm` (full width at half maximum, nm).
#' @keywords datasets
"fwhm.prisma"

#' PRISMA measured spectral response functions
#'
#' Full measured per-band spectral response function (SRF) curves for the
#' PRISMA hyperspectral imager, sampled on a common wavelength grid.
#'
#' @format A data frame with 2091 rows (wavelength grid) and 235 columns:
#'   `wavelength` followed by one column per PRISMA band giving that band's
#'   measured relative spectral response at each wavelength.
#' @keywords datasets
"srf.prisma"

#' Sentinel-2A/2B MSI measured spectral response functions
#'
#' Full measured per-band spectral response function (SRF) curves for the
#' Sentinel-2A and Sentinel-2B MultiSpectral Instrument (MSI), sampled on a
#' common wavelength grid, as distributed by ESA.
#'
#' @format A data frame with 2301 rows (wavelength grid, `SR_WL`) and 14
#'   columns: `SR_WL` followed by one column per MSI band (B1-B12, B8A)
#'   giving that band's measured relative spectral response at each
#'   wavelength.
#' @keywords datasets
#' @name srf-sentinel2
"srf.sentinel2a"

#' @rdname srf-sentinel2
"srf.sentinel2b"

#' Multi-sensor spectral band centre/bounds lookup table
#'
#' A combined lookup table of nominal band bounds and average centre
#' wavelength across the sensors supported by the package's sensor
#' convolution and index functions.
#'
#' @format A data frame with 343 rows (one per sensor/band combination) and
#'   5 columns: `Sensor`, `channel` (band name), `lb`/`ub` (lower/upper
#'   wavelength bound of the band, nm) and `average` (band centre
#'   wavelength, nm).
#' @keywords datasets
"sensor.characteristics"

#' ASTM E490 extraterrestrial solar irradiance spectrum
#'
#' Reference top-of-atmosphere solar spectral irradiance, used to convert
#' simulated top-of-canopy radiance to reflectance and in SPART's
#' atmosphere-coupled simulations.
#'
#' @format A data frame with 2001 rows (400-2400 nm) and 2 columns: `wave`
#'   (wavelength, nm) and `EIrrad` (extraterrestrial spectral irradiance,
#'   W m^-2 nm^-1).
#' @keywords datasets
"Extraterrestrial_irradiance"

#' Angers leaf biochemistry/optics reference database
#'
#' A bundled subset of the Angers leaf-optical-properties experimental
#' database (leaf biochemical traits paired with measured reflectance and
#' transmittance spectra), used as example/validation data for leaf model
#' inversion.
#'
#' @format A named list with components `DataBioch` (data frame of measured
#'   leaf biochemical traits, one row per sample), `lambda` (wavelength
#'   vector for the spectra), `Refl`/`Tran` (measured reflectance and
#'   transmittance spectra, one row per sample) and `nbSamples` (sample
#'   count).
#' @keywords datasets
"LeafDB.Angers"
