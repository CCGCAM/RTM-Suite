Sensor Simulation
=====================

Every model on :doc:`models` simulates a spectrum at native resolution
(400-2500nm, 1nm steps). A real sensor never sees that -- it integrates
incoming light over each band's own spectral response function (SRF):
some wavelengths inside the band contribute fully, others (near the
edges) only partially, and everything outside the band contributes
nothing. **Spectral convolution** is that integration, applied to a
simulated spectrum so it can be compared apples-to-apples with real
sensor data (or used to train a model that will later run on real sensor
bands, as in :doc:`trait_inversion`).

Two ways to convolve, depending on what's known about the sensor
--------------------------------------------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 25 35 40

   * - Function
     - Uses
     - When to use it
   * - :func:`~toolsrtm.srf.spectral_convolution_srf`
     - A real, measured per-wavelength SRF table
     - Sentinel-2A/B (:func:`~toolsrtm.srf.srf_sentinel2a`/
       :func:`~toolsrtm.srf.srf_sentinel2b`), PRISMA
       (:func:`~toolsrtm.srf.srf_prisma`) -- the most accurate option,
       whenever a real measured SRF is available.
   * - :func:`~toolsrtm.srf.spectral_convolution_gaussian`
     - Only nominal band center (+ FWHM or band edges), approximated as a
       Gaussian response
     - EnMAP (:func:`~toolsrtm.srf.enmap_characteristics`), Landsat,
       MODIS, and 9 other bundled sensors
       (:func:`~toolsrtm.srf.sensor_characteristics`), or **any custom
       sensor/camera** -- pass your own ``centers``/``fwhm`` directly.

Sentinel-2A / PRISMA: real measured SRF
--------------------------------------------

.. code-block:: python

   import numpy as np
   from toolsrtm import prospect_d, foursail
   from toolsrtm.srf import srf_sentinel2a, srf_prisma, spectral_convolution_srf

   inputLUT = dict(N=1.5, Cab=40, Car=8, Anth=1, Cbrown=0, EWT=0.01, LMA=0.009, alpha=40,
                    LIDFa=-0.35, LIDFb=-0.15, TypeLidf=1, LAI=3, hspot=0.01, tts=30, tto=0, psi=0)
   sail = foursail(inputLUT, np.full(2101, 0.15), leaf_model="PROSPECT-D", spectrum_all=True)
   wl = np.arange(400, 2501)

   s2a = srf_sentinel2a()
   conv_s2 = spectral_convolution_srf(wl, sail.rsot, s2a)
   print(s2a.band_names, conv_s2.rfl.round(4))
   # ['B1', 'B2', ..., 'B12'] [0.019  0.0247 0.0511 ... ]

   prisma = srf_prisma()
   conv_pr = spectral_convolution_srf(wl, sail.rsot, prisma)
   print(len(conv_pr.rfl), "PRISMA bands")   # 234

Hyperspectral -> multispectral: EnMAP, Landsat, MODIS, custom
--------------------------------------------------------------------

.. code-block:: python

   from toolsrtm.srf import spectral_convolution_gaussian

   conv_enmap = spectral_convolution_gaussian(wl, sail.rsot, sensor="EnMAP")
   print(len(conv_enmap.rfl), "EnMAP channels")   # 242

   conv_modis = spectral_convolution_gaussian(wl, sail.rsot, sensor="MODIS")
   print(len(conv_modis.rfl), "MODIS bands")   # 19

   # Your own sensor/camera: just pass band centers (nm) -- fwhm is optional,
   # approximated from neighbouring-band spacing if omitted.
   my_centers = np.array([550.0, 660.0, 790.0, 1050.0])
   conv_custom = spectral_convolution_gaussian(wl, sail.rsot, centers=my_centers)
   print(conv_custom.rfl.round(4))   # [0.0391 0.0623 0.27   0.3267]

``sensor_characteristics()``'s docstring lists all 13 other bundled
sensors (ALI, Hyperion, Landsat4-8, MODIS, Quickbird, RapidEye,
Sentinel2a/b, WorldView2-4/2-8).

What's next
-----------------

- :doc:`spectral_indices` -- computing vegetation indices from the
  convolved bands above.
- :doc:`trait_inversion` -- training an inversion model on convolved
  bands, so it can later run on real sensor data.
- :doc:`earth_observation` -- applying that model to a real Sentinel-2
  scene via STAC.
