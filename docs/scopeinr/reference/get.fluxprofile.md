# get.fluxprofile

`get.fluxprofile` propagates the top-of-canopy direct and diffuse
irradiance down through the `nl` canopy layers (and reflects them back
up), using the layer transmittance/reflectance factors from
[`get.reflectances`](get.reflectances.md), to obtain the vertical
profile of direct and diffuse radiation fluxes within the canopy.

## Usage

``` r
get.fluxprofile(
  Esun_,
  Esky_,
  rsoil,
  Xss,
  Xsd,
  Xdd,
  R_sd,
  R_dd,
  nl,
  nwl,
  rs.thermal = 0.06
)
```

## Arguments

- Esun\_:

  numeric vector of length `nwl`. Top-of-canopy direct solar irradiance
  spectrum.

- Esky\_:

  numeric vector of length `nwl`. Top-of-canopy diffuse sky irradiance
  spectrum.

- rsoil:

  numeric vector of length `nwl`. Soil reflectance spectrum (bottom
  boundary condition).

- Xss:

  numeric vector of length `nl`. Direct-direct transmittance of each
  layer, as returned by [`get.reflectances`](get.reflectances.md).

- Xsd:

  numeric matrix with `nl` rows and `nwl` columns. Normalized
  direct-to-diffuse transmittance of each layer, as returned by
  [`get.reflectances`](get.reflectances.md).

- Xdd:

  numeric matrix with `nl` rows and `nwl` columns. Normalized
  diffuse-to-diffuse transmittance of each layer, as returned by
  [`get.reflectances`](get.reflectances.md).

- R_sd:

  numeric matrix with `nl+1` rows and `nwl` columns.
  Directional-hemispherical reflectance at the top of each layer, as
  returned by [`get.reflectances`](get.reflectances.md).

- R_dd:

  numeric matrix with `nl+1` rows and `nwl` columns.
  Hemispherical-hemispherical reflectance at the top of each layer, as
  returned by [`get.reflectances`](get.reflectances.md).

- nl:

  integer. Number of canopy layers.

- nwl:

  integer. Number of wavelengths in the spectral domain.

- rs.thermal:

  numeric. Soil reflectance used to initialize the flux matrices and
  (when `Xsd` does not cover the thermal region) to extend the spectral
  domain into the thermal region; default 0.06.

## Value

A list with:

- Es\_:

  numeric matrix with `nl+1` rows and `nwl` columns. Direct solar
  irradiance at the top of each layer (down to the soil).

- Emin\_:

  numeric matrix with `nl+1` rows and `nwl` columns. Downward diffuse
  irradiance at the top of each layer.

- Eplu\_:

  numeric matrix with `nl+1` rows and `nwl` columns. Upward diffuse
  irradiance at the top of each layer.

## Author

     Wout Verhoef, Christiaan van der Tol, Joris Timmermans (Original version in Matlab)

Carlos Camino (Ported version into R)

## Examples

``` r
if (FALSE) { # \dontrun{
Eflux <- get.fluxprofile(Esun_, Esky_, rsoil, Xss, Xsd, Xdd, R_sd, R_dd, nl, nwl)
} # }
```
