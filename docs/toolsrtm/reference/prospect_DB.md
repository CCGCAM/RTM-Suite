# PROSPECT-Dynamic with brown pigments

PROSPECT-Dynamic with brown pigments

## Usage

``` r
prospect_DB(N, Cab, Car, Anth, Cbrown, EWT, LMA, alpha)
```

## Arguments

- N:

  numeric. Leaf structure parameter

- Cab:

  numeric. Chlorophyll content (microg.cm-2)

- Car:

  numeric. Carotenoid content (microg.cm-2)

- Anth:

  numeric. Anthocyanins content (microg.cm-2)

- Cbrown:

  numeric. Brown pigment content (Arbitrary units)

- EWT:

  numeric. Equivalent Water Thickness (g.cm-2). Default is 0.009 Default
  is 0.012

- LMA:

  numeric. Leaf Mass per Area (g.cm-2). Default is 0.012

- alpha:

  numeric. Maximum incidence angle defining the solid angle of incident
  light.By default is 40

## Value

List of lambda with leaf directional-hemisphrical reflectance and
transmittance

## References

Féret J-B, Gitelson AA, Noble SD & Jacquemoud S, 2017. PROSPECT-D:
Towards modeling leaf optical properties through a complete lifecycle.
Remote Sensing of Environment, 193, 204–215.
https://doi.org/10.1016/j.rse.2017.03.004

Jacquemoud S, Baret F, Hanocq J-F, 1992. Modeling spectral and
bidirectional soil reflectance. Remote Sensing of Environment, 41,
123–132. https://doi.org/10.1016/0034-4257(92)90072-R

Jacquemoud, S., Baret, F., 1990. PROSPECT: a model of leaf optical
properties spectra. Remote Sens. Environ. 34, 75–91.
https://doi.org/10.1016/0034-4257 (90)90100-Z.

Authors:

Jean-Baptiste Feret (jb.feret@teledetection.fr)

Frédéric Baret (baret@avignon.inra.fr)

Stephane JAacquemoud (jacquemoud@ipgp.fr)

This function includes numerical optimizations proposed in the FLUSPECT
code

Authors:

Wout Verhoef

Christiaan van der Tol (c.vandertol@utwente.nl)

Joris Timmermans
