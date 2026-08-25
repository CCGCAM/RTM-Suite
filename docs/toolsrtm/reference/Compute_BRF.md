# Computes bidirectional reflectance factor based on outputs from PRO-4SAIL and sun position

Authors of the version:Jean-Baptiste FERET (jb.feret@teledetection.fr)
Florian de BOISSIEU (fdeboiss@gmail.com) Copyright 2019/11 Jean-Baptiste
FERET

## Usage

``` r
Compute_BRF(
  rdot = NULL,
  rsot = NULL,
  tts = NULL,
  data.light = NULL,
  short.waves = T
)
```

## Arguments

- rdot:

  numeric. Hemispherical-directional reflectance factor in viewing
  direction

- rsot:

  numeric. Bi-directional reflectance factor

- tts:

  numeric. Solar zenith angle

- data.light:

  list. direct and diffuse radiation for clear conditions, is NULL use
  default values

- short.waves:

  boolean . Is true the outputs is shorted to wave for Fluspect-model.

## Value

BRF numeric. Bidirectional reflectance factor

## Details

The direct and diffuse light are taken into account as proposed by:
Francois et al. (2002) Conversion of 400-1100 nm vegetation albedo
measurements into total shortwave broadband albedo using a canopy
radiative transfer model, Agronomie

Es = direct Ed = diffuse
