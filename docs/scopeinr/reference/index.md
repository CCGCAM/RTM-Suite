# Package index

## All functions

- [`aggreg()`](aggreg.md) : Aggregates MODTRAN data over SCOPE bands by
  averaging

- [`table_Vcmax_`](biochemistry-tables.md)
  [`table_Jmax_`](biochemistry-tables.md) : Biochemical model lookup
  tables for Vcmax25/Jmax temperature response

- [`brdf_angles`](brdf-angle-datasets.md)
  [`brdf_angles2`](brdf-angle-datasets.md)
  [`brdf_angles_no_oversampling`](brdf-angle-datasets.md) : BRDF
  viewing/illumination angle configurations

- [`calczenithangle()`](calczenithangle.md) : Calculates pi/2-the angle
  of the sun with the slope of the surface.

- [`calc_rssrbs()`](calc_rssrbs.md) : Function for calculating the rss
  and rbs values based on the input parameters (SMC, LAI and rbs)

- [`constants`](constants.md) : Physical and physiological constants
  used by the SCOPE model

- [`count_k()`](count_k.md) : count_k

- [`data.opts`](data.opts.md) : Default SCOPE run options

- [`dcum()`](dcum.md) : dcum function

- [`define.bands()`](define.bands.md) : Define spectral regions used by
  SCOPEinR

- [`define.constants()`](define.constants.md) : Define physical
  constants used by SCOPEinR

- [`define_temp_response_biochem()`](define_temp_response_biochem.md) :
  get the variables needed in temp_response_biochem

- [`get.aggregator.ebal()`](get.aggregator.ebal.md) :
  get.aggregator.ebal

- [`get.BallBerry()`](get.BallBerry.md) :

  get.BallBerry `get.BallBerry` get BerryBall value from Berry Model

- [`get.biochemical.MD12()`](get.biochemical.MD12.md) :
  get.biochemical.MD12

- [`get.biochemical()`](get.biochemical.md) : get.biochemical

- [`get.brdf()`](get.brdf.md) : get.brdf

- [`get.calcTOCirr()`](get.calcTOCirr.md) : get.calcTOCirr

- [`get.calc_fluspect_bcar()`](get.calc_fluspect_bcar.md) : Run the
  combined B/Cx Fluspect leaf model from a packed parameter vector

- [`get.Ci.next()`](get.Ci.next.md) :

  get.Ci.next `get.Ci.next` Function to calculate the difference between
  "guessed" Ci (Ci_in) and Ci computed using BB after computing A
  Test-function for iteration (note that it assigns A in the function's
  context.) As with the next section, this code can be read as if the
  function body executed at this point. (if iteration was used). In
  other words, A is assigned at this point in the file (when iterating).

- [`get.computeA()`](get.computeA.md) :

  get.computeA `get.computeA` Compute the net CO2 assimilation rate
  using the Farquhar model Note: even though computeA() is written as a
  separate function, the code is, in fact, executed exactly this point
  in the file (i.e. between the previous if clause and the next section

- [`get.COST_4Fluspect.for.SCOPE()`](get.COST_4Fluspect.for.SCOPE.md) :
  Cost function for fitting Fluspect leaf parameters to measured
  reflectance/transmittance

- [`get.e2phot()`](get.e2phot.md) : get.e2phot calculates the number of
  moles of photons

- [`get.ebal()`](get.ebal.md) :

  get.ebal `get.ebal` Calculates the energy balance of a vegetated
  surface

- [`get.ephoton()`](get.ephoton.md) : get.ephoton

- [`get.Fluorescence.model()`](get.Fluorescence.model.md) :

  get.Fluorescence.model `get.Fluorescence.model` Fluorescence model

- [`get.fluspect_mSCOPE()`](get.fluspect_mSCOPE.md) :
  get.fluspect_mSCOPE

- [`get.fluxprofile()`](get.fluxprofile.md) : get.fluxprofile

- [`get.gsFun()`](get.gsFun.md) :

  get.gsFun `get.gsFun` get stomatal conductance

- [`get.heatfluxes()`](get.heatfluxes.md) :

  get.heatfluxes `get.heatfluxes` Calculates latent and sensible heat
  flux

- [`get.high.temp.inhibtionC3()`](get.high.temp.inhibtionC3.md) :

  get.high.temp.inhibtionC3 `get.high.temp.inhibtionC3` High Temperature
  Inhibition Function:The following function pertains to C3
  photosynthesis

- [`get.MD12()`](get.MD12.md) : MD12 algorithm for the computation of
  fluorescence yield

- [`get.merge.SCOPE()`](get.merge.SCOPE.md) : Get.merge.Output generate
  the LUT table adding the apparent reflectance, radiance or
  fluorescence emission

- [`get.Monin.Obukhov()`](get.Monin.Obukhov.md) : get.Monin.Obukhov
  function

- [`get.numjacobian()`](get.numjacobian.md) : Compute the numerical
  Jacobian of the Fluspect leaf model with respect to its parameters

- [`get.outs.in()`](get.outs.in.md) :

  `get.outs.in` get the outputs from the SCOPE model. This code is for
  single simulations

- [`get.outs.lut()`](get.outs.lut.md) :

  `get.outs.lut` get the outputs from SCOPE model. This code is for a
  LUT of simulations

- [`get.outs.lut.v2()`](get.outs.lut.v2.md) :

  `get.outs.lut` get the outputs from SCOPE model. This code is for a
  LUT of simulations

- [`get.phstar()`](get.phstar.md) : subfunction phs for stability
  correction (eg. Paulson, 1970)

- [`get.Planck()`](get.Planck.md) : Planck function

- [`get.psih()`](get.psih.md) : subfunction ph for stability correction
  (eg. Paulson, 1970)

- [`get.psim()`](get.psim.md) : subfunction pm for stability correction
  (eg. Paulson, 1970)

- [`get.Pso()`](get.Pso.md) : get.Pso function

- [`get.reflectances()`](get.reflectances.md) : get.reflectances

- [`get.resistances()`](get.resistances.md) : get.resistances

- [`get.RTMf()`](get.RTMf.md) :

  get.RTMf `get.RTMf` Calculates the spectrum of fluorescent radiance in
  the observer's direction and also the TOC spectral hemispherical
  upward Fs flux.

- [`get.RTMt.planck()`](get.RTMt.planck.md) :

  get.RTMt.planck `get.RTMt.planck`analogue to get.RTMt.sb, this
  function calculates total outgoing radiation in hemispherical
  direction and total absorbed radiation per leaf and soil component.
  Radiation is integrated over the whole thermal spectrum with
  Stefan-Boltzman's equation. This function is a simplified version of
  'get.RTMt.planck', and is less time consuming since it does not do the
  calculation for each wavelength separately.

- [`get.RTMt.sb()`](get.RTMt.sb.md) :

  getRTMt.sb `get.RTMt.sb` Calculates total outgoing radiation in
  hemispherical direction and total absorbed radiation per leaf and soil
  component. Radiation is integrated over the whole thermal spectrum
  with Stefan-Boltzman's equation. This function is a simplified version
  of 'getRTMt.planck', and is less time consuming since it does not do
  the calculation for each wavelength separately.

- [`get.RTMz()`](get.RTMz.md) :

  get.RTMz `get.RTMz` Calculates the small modification of TOC outgoing
  radiance due to the conversion of Violaxanthin into Zeaxanthin in
  leaves

- [`get.SCOPE.ind()`](get.SCOPE.ind.md) : get simulations based on SCOPE
  model

- [`get.SCOPE.outputs()`](get.SCOPE.outputs.md) :

  `get.SCOPE.outputs` get the outputs from SCOPE model

- [`get.SCOPE.parallel()`](get.SCOPE.parallel.md) : Run SCOPE
  simulations in parallell

- [`get.SCOPE.plots()`](get.SCOPE.plots.md) :

  `get.SCOPE.plots` get simulations based on SCOPE model

- [`get.SCOPE()`](get.SCOPE.md) :

  `get.SCOPE` get simulations based on SCOPE model

- [`get.spectra.SCOPE()`](get.spectra.SCOPE.md) :

  get.spectra.SCOPE `get.spectra.SCOPE` Calculates the spectra of
  hemisperical and directional observed \#' a function to get the
  spectral characteristics for SCOPE model

- [`get.Stefan_Boltzmann()`](get.Stefan_Boltzmann.md) : Stefan-Boltzmann
  equation

- [`get.temperature.functionC3()`](get.temperature.functionC3.md) :

  get.temperature.functionC3 `get.temperature.functionC3` Temperature
  Correction Functions:The following function pertains to C3
  photosynthesis

- [`get.volscatt.scope()`](get.volscatt.scope.md) : get.volscatt.scope
  version 2.0 from SCOPE model

- [`get.zo_and_d()`](get.zo_and_d.md) :

  get.zo_and_d model `get.zo_and_d` Calculates roughness length for
  momentum and zero plane displacement from vegetation height and LAI

- [`getBSM()`](getBSM.md) : Brightness-Shape-Moisture soil model

- [`getCSV()`](getCSV.md) : Get CSV Data from Multiple Folders

- [`getFluspect.B.SCOPE()`](getFluspect.B.SCOPE.md) : Leaf FLUSPECT-B
  model for SCOPE

- [`getFluspect.Cx.SCOPE()`](getFluspect.Cx.SCOPE.md) : Leaf
  FLUSPECT-B-Cx model for SCOPE

- [`getinputLUT()`](getinputLUT.md) : Get main input for SCOPE

- [`getLUT.SCOPE()`](getLUT.SCOPE.md) : Generate LUT for SCOPE

- [`getLUT.SCOPE.v1()`](getLUT.SCOPE.v1.md) : Generate LUT for SCOPE
  (v1)

- [`getLUT_time()`](getLUT_time.md) : Get LUT table for SCOPE model

- [`getRTMo()`](getRTMo.md) :

  `getRTMo` Calculates the spectra of hemisperical and directional
  observed visible and thermal radiation (fluxes E and radiances L), as
  well as the single and bi-directional gap probabilities

- [`latin_hypercube_input()`](latin_hypercube_input.md) :
  latin_hypercube_input function

- [`leafangles()`](leafangles.md) : Subroutine FluorSail_dladgen
  (Version 2.3)

- [`leaf_spectrum`](leaf_spectrum.md) : Example leaf
  reflectance/transmittance spectrum

- [`MD12()`](MD12.md) : MD12 algorithm for the computation of
  fluorescence yield

- [`meanleaf()`](meanleaf.md) : Calculates the layer average and the
  canopy average of leaf properties per layer, per leaf angle and per
  leaf azimuth (36)

- [`meanleaf.v2()`](meanleaf.v2.md) : meanleaf.v2

- [`Rin_`](meteo-timeseries.md) [`Rli_`](meteo-timeseries.md)
  [`Esun_`](meteo-timeseries.md) [`Esky_`](meteo-timeseries.md)
  [`Ta_`](meteo-timeseries.md) [`ea_`](meteo-timeseries.md)
  [`p_`](meteo-timeseries.md) [`u_`](meteo-timeseries.md)
  [`t_`](meteo-timeseries.md) [`year_`](meteo-timeseries.md) : Example
  meteorological/radiation time series for diurnal SCOPE simulations

- [`optipar`](optipar-datasets.md) [`optipar.2015`](optipar-datasets.md)
  [`optipar2017.ProspectD`](optipar-datasets.md)
  [`optipar2020.prospectD.BSM2019`](optipar-datasets.md)
  [`optipar2021.Pro.CX`](optipar-datasets.md) : Leaf optical parameters
  (PROSPECT/Fluspect, various parameterizations)

- [`satvap()`](satvap.md) : calculates the saturated vapour pressure at
  temperature T (degrees C) and the derivative of es to temperature s
  (kPa/C)

- [`inputsSCOPE`](scope-input-datasets.md)
  [`input_border`](scope-input-datasets.md) : Example SCOPE input LUT
  and input border/range definitions

- [`SCOPE.LUT.default`](SCOPE.LUT.default.md) : Default SCOPE input
  look-up table

- [`sel_root()`](sel_root.md) : quadratic formula, root of least
  magnitude

- [`Sint()`](Sint.md) : Simpson-like trapezoidal integration

- [`slope_satvap()`](slope_satvap.md) : calculates the saturated vapour
  pressure at temperature T (degrees C) and the derivative of es to
  temperature s (kPa/C)

- [`soil.rfl`](soil.rfl.md) : Example soil reflectance spectrum

- [`soilwat()`](soilwat.md) :

  soilwat function `soilwat` In this model it is assumed that the water
  film area is built up

- [`Soil_Inertia0()`](Soil_Inertia0.md) : Calculate the soil thermal
  inertia from known soil thermal properties

- [`Soil_Inertia1()`](Soil_Inertia1.md) : Soil thermal inertia method by
  Murray and Verhoef

- [`soil_respiration()`](soil_respiration.md) : soil respiration

- [`tav()`](tav.md) : Stern's formula in Lekner & Dorf (1988) gives
  reflectance for alfa = 90 degrees
