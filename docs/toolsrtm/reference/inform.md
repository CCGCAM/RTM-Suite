# The INvertible FOrest Reflectance Model coupled with leaf models

`inform` INFORM simulation based on a set of combinations of
inpuparameters

## Usage

``` r
inform(inputLUT = NULL, rsoil = rsoil, LeafModel = "PROSPECT-PRO")
```

## Arguments

- inputLUT:

  LUT table with distribution of biophysical parameters used as input
  parameters in the model

- rsoil:

  numeric. Soil reflectance

- LeafModel:

  Version of PROSPECT model; 'PRO' or 'D' is accepted. By default
  'PROSPECT-PRO' is used. Liberty model and Fluspect-B and Fluspect-B-Cx
  is also provided. Options are : 'PROSPECT-PRO','PROSPECT-D',
  'Liberty', 'Fluspect-B' and 'Fluspect-B-Cx'.

## Value

rfl

## References

INFORM (Atzberger, 2000; Schlerf&Atzberger, 2006) simulates the
bi-directional reflectance of forest stands between 400 and 2500 nm.

INFORM is essentially an innovative combination of FLIM (Rosema et al.,
1992), SAIL (Verhoef, 1984), and PROSPECT (Jacquemoud et al. 1996)

Atzberger, C. 2000: Development of an invertible forest reflectance
model: The INFOR-Model. In: Buchroithner (Ed.): A decade of
trans-european remote sensing cooperation. Proceedings of the 20th
EARSeL Symposium Dresden, Germany, 14.-16. June 2000: 39-44.

Schlerf, M. & Atzberger, C. (2006): Inversion of a forest reflectance
model to estimate biophysical canopy variables from hyperspectral remote
sensing data. Remote Sensing of Environment, 100: 281-294

Rosema, A., Verhoef, W., Noorbergen, H. 1992: A new forest light
interaction model in support of forest monitoring. Remote Sensing of
Environment, 42: 23-41.

Jacquemoud S., Ustin S.L., Verdebout J., Schmuck G., Andreoli G.,
Hosgood B. (1996): Estimating leaf biochemistry using the PROSPECT leaf
optical properties model, Remote Sens. Environ., 56:194-202.

Verhoef, W. 1984: Light scattering by leaf layers with application to
canopy reflectance modeling: The SAIL model. Remote Sensing of
Environment, 16: 125-141.

Basic version of INFORM: Clement Atzberger, 1999 INFORM modifications
and validation: Martin Schlerf, 2004-2007

## Author

Clement Atzberger and Martin Schlerf (Original version in Matlab and
python)

Carlos Camino (Ported version into R)
