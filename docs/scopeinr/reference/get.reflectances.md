# get.reflectances

`get.reflectances` propagates the thin-layer reflectance and
transmittance factors of the canopy (computed by the SAIL 4-stream
radiative transfer scheme) down through all `nl` canopy layers to the
soil, and back up, to obtain the directional-hemispherical and
hemispherical-hemispherical reflectance at the top of each layer.

## Usage

``` r
get.reflectances(tau_ss, tau_sd, tau_dd, rho_dd, rho_sd, rsoil, nl, nwl = 2162)
```

## Arguments

- tau_ss:

  numeric matrix with `nl` rows and `nwl` columns. Direct-direct
  transmittance of each thin canopy layer.

- tau_sd:

  numeric matrix with `nl` rows and `nwl` columns. Direct-diffuse
  transmittance of each thin canopy layer.

- tau_dd:

  numeric matrix with `nl` rows and `nwl` columns. Diffuse-diffuse
  transmittance of each thin canopy layer.

- rho_dd:

  numeric matrix with `nl` rows and `nwl` columns. Diffuse-diffuse
  reflectance of each thin canopy layer.

- rho_sd:

  numeric matrix with `nl` rows and `nwl` columns. Direct-diffuse
  reflectance of each thin canopy layer.

- rsoil:

  numeric vector of length `nwl`. Soil reflectance spectrum at the
  bottom boundary.

- nl:

  integer. Number of canopy layers.

- nwl:

  integer. Number of wavelengths in the spectral domain (default 2162,
  i.e. 400-2500 nm plus the thermal region).

## Value

A list with:

- R_sd:

  numeric matrix with `nl+1` rows and `nwl` columns.
  Directional-hemispherical reflectance at the top of each layer
  (including the soil, layer nl+1).

- R_dd:

  numeric matrix with `nl+1` rows and `nwl` columns.
  Hemispherical-hemispherical reflectance at the top of each layer
  (including the soil, layer nl+1).

- Xss:

  numeric vector of length `nl`. Direct-direct transmittance of each
  layer (identical to `tau_ss`, kept for use in flux propagation).

- Xsd:

  numeric matrix with `nl` rows and `nwl` columns. Normalized
  direct-to-diffuse transmittance of each layer, accounting for multiple
  reflections with the layers below.

- Xdd:

  numeric matrix with `nl` rows and `nwl` columns. Normalized
  diffuse-to-diffuse transmittance of each layer, accounting for
  multiple reflections with the layers below.

## Author

     Wout Verhoef, Christiaan van der Tol, Joris Timmermans (Original version in Matlab)

Carlos Camino (Ported version into R)

## Examples

``` r
if (FALSE) { # \dontrun{
refl <- get.reflectances(tau_ss, tau_sd, tau_dd, rho_dd, rho_sd, rsoil, nl, nwl = 2162)
} # }
```
