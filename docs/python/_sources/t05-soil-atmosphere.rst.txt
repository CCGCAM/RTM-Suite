05. Soil & Atmosphere
==========================

What you will learn
------------------------

- What MARMIT, BSM, SMAC, and SPART each simulate.
- How soil moisture changes a soil reflectance spectrum, and how much
  that matters relative to the canopy above it.
- Why top-of-canopy (TOC) and top-of-atmosphere (TOA) reflectance
  diverge sharply at specific wavelengths.

Concept
-----------

Two independent soil models feed two different canopy/atmosphere
chains:

.. list-table::
   :header-rows: 1
   :widths: 20 45 35

   * - Model
     - What it represents
     - Used by
   * - MARMIT
     - Starts from a real *dry* reference soil spectrum and adds a
       physically modelled liquid-water film, so the same soil can be
       simulated at any moisture level.
     - Standalone, or as ``rsoil`` into any canopy model
       (:doc:`t04-canopy-models`).
   * - BSM (Brightness-Shape-Moisture)
     - Builds a soil spectrum from three empirical parameters
       (``BSMBrightness``, ``BSMlat``, ``BSMlon``) plus a wetting term --
       no reference spectrum needed.
     - Built into :func:`~toolsrtm.spart.spart_toa`.
   * - SMAC
     - Atmospheric radiative transfer (gas absorption + aerosol
       scattering) converting TOC reflectance into what a sensor
       measures above the atmosphere.
     - Built into :func:`~toolsrtm.spart.spart_toa`.
   * - SPART
     - The full chain: BSM soil -> fourSAIL canopy -> SMAC atmosphere ->
       TOA reflectance, already resampled to a real sensor's bands.
     - :func:`~toolsrtm.spart.spart_toa`.

Python tools used
----------------------

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Function
     - Key arguments
   * - :func:`~toolsrtm.marmit.get_marmit_rsoil`
     - ``soil_id`` (dry reference spectrum), ``L`` (water-film thickness,
       cm), ``eps`` (roughness). Returns ``.rsoil_dry``/``.rsoil_wet``
       (2101-element each) and ``.smc`` (derived soil moisture content).
   * - :func:`~toolsrtm.spart.spart_toa`
     - Trait dict (leaf + canopy + atmosphere: ``Pa, aot550, uo3, uh2o``),
       ``sensor`` (a :class:`~toolsrtm.smac.SmacSensor`, e.g.
       :func:`~toolsrtm.smac.sentinel2a_msi`), ``leaf_model``, plus BSM
       soil kwargs (``BSMBrightness, BSMlat, BSMlon, SMp``). Returns
       ``.wl_smac`` (band centers), ``.rfl_toc_brdf``, ``.rfl_toa`` --
       one value per sensor band, not the full 1nm spectrum.

Run the example
--------------------

.. code-block:: python

   from toolsrtm import get_marmit_rsoil, spart_toa, sentinel2a_msi

   # 1. MARMIT: dry -> wet soil spectrum
   soil = get_marmit_rsoil(soil_id=3, L=0.05, eps=0.4, version="marmit1")
   print("Soil moisture content:", float(soil.smc))
   for wl in (550, 850, 1600):
       i = wl - 400
       print(f"{wl}nm: dry={soil.rsoil_dry[i]:.4f}  wet={soil.rsoil_wet[i]:.4f}")

   # 2. SPART: BSM soil -> fourSAIL canopy -> SMAC atmosphere -> TOA
   result = spart_toa(
       dict(N=1.5, Cab=40, Car=8, Anth=1, Cbrown=0, EWT=0.01, LMA=0.009, alpha=40,
            Prot=0.002, CBC=0.007, LIDFa=-0.35, LIDFb=-0.15, TypeLidf=1,
            LAI=3, hspot=0.01, tts=30, tto=0, psi=0,
            Pa=1000, aot550=0.3246, uo3=0.348, uh2o=1.4116),
       sensor=sentinel2a_msi(), leaf_model="PROSPECT-PRO",
       BSMBrightness=0.5, BSMlat=25, BSMlon=45, SMp=15,
   )
   print("Sentinel-2A bands (nm):", result.wl_smac)
   print("TOC:", result.rfl_toc_brdf.round(4))
   print("TOA:", result.rfl_toa.round(4))

Result
----------

Printed output (exact, deterministic)::

   Soil moisture content: 39.900691712115204
   550nm: dry=0.1316  wet=0.1064
   850nm: dry=0.2608  wet=0.2151
   1600nm: dry=0.4758  wet=0.3178
   Sentinel-2A bands (nm): [ 445.  520.  560.  654.  701.  743.  779.  789.  871.  942. 1372. 1639. 2256.]
   TOC: [0.0159 0.0414 0.0556 0.0238 0.0662 0.3305 0.4021 0.4032 0.4064 0.4029 0.2806 0.2248 0.0805]
   TOA: [0.1182 0.1057 0.0899 0.0488 0.0803 0.3104 0.3846 0.3679 0.3904 0.1215 0.0022 0.2091 0.0736]

.. figure:: _figures/marmit_soil.png
   :alt: Dry vs wet soil reflectance spectrum from MARMIT, real output of the code above
   :width: 75%

   Real output: MARMIT's dry-reference vs. wetted soil reflectance -- the SWIR water-absorption dips (~1400/1900nm) deepen and overall brightness drops as the soil wets.

.. figure:: _figures/spart_toc_toa.png
   :alt: SPART TOC and TOA reflectance across Sentinel-2A bands, real output of the code above
   :width: 75%

   Real output: TOC and TOA reflectance diverge sharply in the water-vapour bands (~940/1370nm), where the atmosphere absorbs most of the signal before it reaches the sensor.

Interpretation
-------------------

Wetting the soil (MARMIT) lowers reflectance everywhere, but not by an
equal fraction: at 1600nm (a SWIR water-absorption region) the drop is
large (0.476 -> 0.318, -33%), while at 550nm (visible) it's smaller
(0.132 -> 0.106, -19%) -- soil moisture leaves its strongest fingerprint
in the SWIR, the same spectral region leaf equivalent-water-thickness
(``EWT``, :doc:`t02-parameters-traits`) affects, for the same physical
reason (liquid water absorbs strongly there).

The TOC-vs-TOA comparison from SPART tells a different, atmosphere-driven
story: at most bands TOC and TOA are reasonably close (e.g. 779nm: 0.403
vs. 0.385), but at the two bands SPART itself sits on real water-vapour
absorption features (942nm, 1372nm), the atmosphere removes most of the
signal before it reaches "space" -- TOC 0.403 collapses to TOA 0.122 at
942nm, and TOC 0.281 collapses to a near-zero 0.002 at 1372nm. This is
exactly why real sensors either avoid placing bands there (Sentinel-2's
narrow B9/B10 bands specifically target these features for atmospheric
correction, not surface retrieval) or need real atmospheric correction
before the data is usable for vegetation analysis.

Try it yourself
--------------------

- Raise ``L`` (MARMIT) from 0.05 to 0.15 (much wetter) and see how much
  further the SWIR dips deepen.
- Change ``aot550`` (aerosol optical thickness) from 0.32 to 0.05 (clear
  sky) and compare how much the visible-band TOA values shift -- aerosol
  scattering matters most in the blue/visible, not the SWIR.
- Feed MARMIT's ``soil.rsoil_wet`` directly as ``rsoil`` into
  :doc:`t04-canopy-models`'s ``foursail()`` call, and compare the
  resulting TOC reflectance against the flat ``rsoil=0.15`` baseline
  used there.

Common mistakes
--------------------

- MARMIT and BSM are two *independent* soil models -- ``spart_toa``
  always uses BSM internally (``BSMBrightness``/``BSMlat``/``BSMlon``/
  ``SMp``), not whatever MARMIT spectrum you may have built separately.
- ``spart_toa`` returns one value **per sensor band** (13 for
  Sentinel-2A), not the full 2101-element native spectrum --
  ``rfl_toc_brdf``/``rfl_toa`` are already resampled.
- Atmospheric water-vapour bands can swing from a plausible reflectance
  value at TOC to near-zero at TOA -- a real atmospheric effect, not a
  sign the simulation broke.

Next
--------

:doc:`t06-scope` -- SCOPE, the ecosystem-scale model that replaces the
"soil brightness + fixed temperature" shortcuts above with a real
coupled energy balance.

----

Using R? -> `ToolsRTM Tutorial 03: SPART
<https://ccgcam.github.io/RTM-Suite/toolsrtm/articles/t03-spart.html>`_
and `Tutorial 16: MARMIT + fourSAIL + SPART
<https://ccgcam.github.io/RTM-Suite/toolsrtm/articles/t16-marmit-soil-in-canopy.html>`_
