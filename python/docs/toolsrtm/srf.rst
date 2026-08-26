toolsrtm.srf
==============

Sensor convolution for the cases ``smac`` doesn't cover: a plain per-band
spectral response function (SRF) with no atmospheric-correction coefficients
(PRISMA, Sentinel-2A/B), and nominal center+FWHM-only convolution (EnMAP,
Landsat, MODIS, and any custom sensor or camera). Direct port of
``ToolsRTM::get.spectral.convolution.srf``/``get.spectral.convolution.gaussian``.

Quick example
-------------

.. code-block:: python

   import numpy as np
   from toolsrtm import foursail
   from toolsrtm.srf import srf_sentinel2a, spectral_convolution_srf

   inputLUT = dict(N=1.5, Cab=40, Car=8, Anth=1, Cbrown=0, EWT=0.01, LMA=0.009, alpha=40,
                    LIDFa=-0.35, LIDFb=-0.15, TypeLidf=1,
                    LAI=3, hspot=0.01, tts=30, tto=0, psi=0)
   sail = foursail(inputLUT, np.full(2101, 0.15), leaf_model="PROSPECT-D", spectrum_all=True)
   wl = np.arange(400, 2501)

   s2a = srf_sentinel2a()
   conv = spectral_convolution_srf(wl, sail.rsot, s2a)
   print(s2a.band_names, conv.rfl.round(4))

.. code-block:: text

   Input                              spectral_convolution_srf() Output
   ---------------------------        ----------------------     ---------------------------
   wl   [n]  native wavelengths                                  conv.wl    band-center wl [b]
   refl [n]  native 1nm spectrum      -------------------->       conv.rfl   convolved refl. [b]
   sensor = srf_sentinel2a() (or a                                conv.band_names
            custom SRF/Gaussian sensor)

.. automodule:: toolsrtm.srf
   :members:
   :undoc-members:
   :show-inheritance:
