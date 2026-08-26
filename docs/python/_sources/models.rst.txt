Radiative Transfer Models
============================

Every function on this site belongs to one of a small number of named
physical models. This page is what each one *is* -- scale, what it
simulates, and how it differs from its siblings -- before :doc:`workflows`
shows how they chain together and :doc:`examples` runs them.

At a glance
--------------

.. list-table::
   :header-rows: 1
   :widths: 15 12 20 25 28

   * - Model
     - Scale
     - Simulates
     - Typical use
     - Python entry point
   * - PROSPECT-D / PROSPECT-PRO
     - Leaf
     - Reflectance/transmittance (R/T)
     - Pigments, water, dry matter (PRO: + protein, CBC)
     - :func:`~toolsrtm.leaf.prospect_d` / :func:`~toolsrtm.leaf.prospect_pro`
   * - Fluspect-B / Fluspect-Cx
     - Leaf
     - R/T + fluorescence
     - Solar-induced fluorescence (SIF)
     - :func:`~toolsrtm.fluspect.fluspect_b` / :func:`~toolsrtm.fluspect.fluspect_cx`
   * - LIBERTY
     - Leaf (needle)
     - R/T
     - Conifer canopies, not broadleaf
     - :func:`~toolsrtm.liberty.liberty`
   * - fourSAIL
     - Canopy
     - TOC reflectance (BRDF)
     - Crops, grassland, single-layer canopies
     - :func:`~toolsrtm.canopy.foursail`
   * - fourSAIL2
     - Canopy
     - TOC reflectance (BRDF)
     - Mixed green + senescent/brown canopies
     - ``toolsrtm.canopy.foursail2``
   * - INFORM
     - Canopy (forest stand)
     - TOC reflectance (BRDF)
     - Forests -- explicit crown/gap geometry
     - :func:`~toolsrtm.inform.inform`
   * - MARMIT
     - Soil
     - Reflectance
     - Soil moisture, from a real dry reference spectrum
     - :func:`~toolsrtm.marmit.get_marmit_rsoil`
   * - BSM
     - Soil
     - Reflectance
     - Soil moisture, from empirical brightness/shape parameters
     - :func:`~scopeinpython.soil.get_bsm`
   * - SPART
     - Soil -> canopy -> atmosphere
     - TOC + TOA reflectance
     - A vegetated scene as a real sensor would measure it
     - :func:`~toolsrtm.spart.spart_toa`
   * - SCOPE
     - Ecosystem (leaf -> canopy -> atmosphere-facing)
     - Reflectance + SIF + energy balance
     - Physiology: temperature, photosynthesis, carbon/water flux
     - :func:`~scopeinpython.scope.get_scope`

This table is the "which model do I want" lookup; the sections below are
the "what does this model actually do" reference for each row.

Leaf optical models
------------------------

All leaf models take pigment/water/dry-matter traits and return a
reflectance/transmittance spectrum -- see the :doc:`t02-parameters-traits` for exactly
what each trait (``N``, ``Cab``, ``EWT``, ...) means.

.. list-table::
   :header-rows: 1
   :widths: 20 45 35

   * - Model
     - What it represents
     - How it differs from the others
   * - :func:`~toolsrtm.leaf.prospect_d` (PROSPECT-D)
     - The reference leaf model: a stack of ``N`` absorbing/scattering
       plates, pigments (``Cab``/``Car``/``Anth``) + water (``EWT``) +
       one lumped dry-matter term (``LMA``).
     - The default -- broadleaf, no fluorescence, dry matter as a single
       term.
   * - :func:`~toolsrtm.leaf.prospect_pro` (PROSPECT-PRO)
     - Same physics as PROSPECT-D, but splits ``LMA`` into ``Prot``
       (protein) + ``CBC`` (cellulose+lignin).
     - Only differs from PROSPECT-D in how dry matter is parameterized --
       useful when protein/nitrogen matters on its own.
   * - :func:`~toolsrtm.fluspect.fluspect_b` (Fluspect-B)
     - PROSPECT-D's optics *plus* chlorophyll-fluorescence excitation-
       emission matrices (``MbI``/``MbII``) -- needed wherever SIF is
       simulated.
     - Adds fluorescence on top of PROSPECT-D; reflectance/transmittance
       themselves are near-identical to PROSPECT-D.
   * - :func:`~toolsrtm.fluspect.fluspect_cx` (Fluspect-B-Cx)
     - Fluspect-B plus a ``Cx`` (xanthophyll de-epoxidation / NPQ) term,
       letting fluorescence yield respond to photoprotection state, not
       just pigment content.
     - The only leaf model with a photoprotection (``Cx``) term; what
       SCOPE's own leaf-optics step uses internally.
   * - :func:`~toolsrtm.liberty.liberty` (LIBERTY)
     - A structurally different model built for conifer needles (Dawson
       et al. 1998) -- explicit cell diameter/intercellular air space
       instead of PROSPECT's Fresnel-refraction layer.
     - Not a PROSPECT variant at all; needle-specific anatomy, a
       genuinely different internal structure, not just different
       defaults.

Canopy models
------------------

All three take a leaf model's reflectance/transmittance plus a soil
background and canopy structure, and return top-of-canopy (TOC) BRDF.

.. list-table::
   :header-rows: 1
   :widths: 20 45 35

   * - Model
     - What it represents
     - How it differs from the others
   * - :func:`~toolsrtm.canopy.foursail` (fourSAIL)
     - The classic PROSAIL turbid-medium canopy: a single, statistically
       homogeneous "cloud" of leaves at a given LAI and angle
       distribution -- no explicit 3D structure.
     - The reference/default; single-layer, single leaf-biochemistry
       profile.
   * - ``toolsrtm.canopy.foursail2`` (fourSAIL2)
     - Two-layer canopy (a green layer + a brown/senescent layer,
       ``fraction_brown``-weighted) -- e.g. a canopy with visible
       dead/dry material mixed in with live foliage.
     - Same turbid-medium idea as fourSAIL, but two vertically-stacked
       layers instead of one.
   * - :func:`~toolsrtm.inform.inform` (INFORM)
     - Forest-stand extension (Atzberger): explicit tree crowns (stem
       density, crown diameter, height) over an understorey +
       background, rather than one homogeneous canopy.
     - The only one of the three with real gap/shadow geometry --
       produces visibly lower reflectance than fourSAIL at the same LAI,
       matching a discontinuous forest stand's real physics.

Soil
--------

.. list-table::
   :header-rows: 1
   :widths: 20 45 35

   * - Model
     - What it represents
     - How it differs from the others
   * - :func:`~toolsrtm.marmit.get_marmit_rsoil` (MARMIT)
     - Starts from a real *dry* reference soil spectrum and adds a
       physically modelled liquid-water film, so the same soil can be
       simulated at any moisture level.
     - The only soil model driven by an actual measured reference
       spectrum rather than empirical shape parameters.
   * - :func:`~scopeinpython.soil.get_bsm` (BSM, Brightness-Shape-Moisture)
     - Builds a soil spectrum from three empirical parameters
       (``BSMBrightness``, ``BSMlat``, ``BSMlon``) plus a wetting term --
       no reference spectrum needed.
     - Purely parametric, not spectrum-driven -- SPART's and SCOPE's own
       soil model.

Soil-Plant-Atmosphere
--------------------------

.. list-table::
   :header-rows: 1
   :widths: 20 45 35

   * - Model
     - What it represents
     - How it differs from the others
   * - :func:`~toolsrtm.smac.sentinel2a_msi` + ``get_smac`` (SMAC)
     - Atmospheric radiative transfer (gas absorption + aerosol
       scattering) that converts top-of-canopy reflectance into what a
       real satellite sensor would measure above the atmosphere.
     - Not a soil or canopy model -- the atmosphere step, only relevant
       when going all the way to top-of-atmosphere (TOA).
   * - :func:`~toolsrtm.spart.spart_toa` (SPART)
     - Not a new physical model, but the full chain: BSM soil ->
       fourSAIL canopy -> SMAC atmosphere -> TOA reflectance, already
       resampled to a real sensor's bands in one call.
     - The end-to-end soil-plant-atmosphere pipeline; no separate
       "simulate native, then convolve" step, unlike plain fourSAIL.

Energy balance / fluorescence: SCOPE
-----------------------------------------

SCOPE (:func:`scopeinpython.scope.get_scope`) is a different kind of
model from everything above, not just a bigger one: instead of assuming
leaf/soil temperature and computing reflectance alone, it **iteratively
solves** leaf and soil temperature so that absorbed radiation balances
sensible + latent heat + photosynthesis (the energy balance), then
derives fluorescence and carbon flux from that solved state. Five
distinct components, chained together:

.. list-table::
   :header-rows: 1
   :widths: 20 45 35

   * - Component
     - What it does
     - Depends on
   * - :func:`~scopeinpython.fluspect.get_fluspect_cx_scope` (+ :func:`~scopeinpython.fluspect_mscope.fluspect_mscope` for multi-layer)
     - Leaf optics: reflectance/transmittance + fluorescence excitation-
       emission matrices, per canopy layer. Same PROSPECT/Fluspect
       physics as above -- just the SCOPE-specific wrapper.
     - Leaf biochemistry traits (:doc:`t02-parameters-traits`)
   * - :func:`~scopeinpython.rtmo.run_rtmo`
     - Optical top-of-canopy BRDF -- physically the same turbid-medium
       idea as fourSAIL, re-implemented to plug into the layers below.
       Stops at reflectance; assumes no temperature yet.
     - Fluspect-Cx output + canopy structure + BSM soil
   * - :func:`~scopeinpython.ebal.ebal`
     - The energy balance: iterates leaf/soil temperature until the
       flux budget closes, calling the biochemistry model at every
       candidate temperature.
     - RTMo output + aerodynamic resistances + meteorology
   * - :func:`~scopeinpython.biochemical.get_biochemical`
     - Leaf-level photosynthesis (Farquhar/Collatz) and fluorescence
       yield, given a leaf micro-environment.
     - Photosynthesis + NPQ traits (:doc:`t02-parameters-traits`)
   * - :func:`~scopeinpython.rtmf.rtmf` / :func:`~scopeinpython.rtmz.rtmz` (optional)
     - Canopy-level fluorescence radiance/flux, and a small zeaxanthin
       (photoprotection) correction to the TOC spectrum.
     - ebal output + Fluspect-Cx's fluorescence matrices

What's next
-----------------

- :doc:`workflows` -- how these models chain together into the four
  standard simulation pipelines.
- :doc:`t02-parameters-traits` -- what every input trait to these models means, its
  unit, and its realistic range.
- :doc:`examples` -- every model above, run with real, verified code.
