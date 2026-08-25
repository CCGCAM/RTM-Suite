# get.fluspect_mSCOPE

`get.fluspect_mSCOPE` is an adaptation of the Fluspect leaf model for
the multi-layer (mSCOPE) canopy representation: leaf optical properties
(reflectance, transmittance, fluorescence excitation-emission matrices
`Mb`/`Mf`, and pigment contribution factors) are computed once per
distinct leaf biochemistry profile layer (`mly$nly` layers), and then
replicated across the `nl` canopy layers that fall within each profile
layer, weighted by `mly$pLAI`.

## Usage

``` r
get.fluspect_mSCOPE(
  mly,
  spectral,
  leafbio,
  soil,
  optipar,
  nl,
  step,
  get.plots = T
)
```

## Arguments

- mly:

  list. Multi-layer leaf biochemistry profile, with elements `nly`
  (number of distinct biochemistry layers), `pLAI` (LAI fraction of each
  layer), and per-layer biochemistry vectors `pCab`, `pEWT`, `pCar`,
  `pLMA`, `pCs`, `pN` (one value per layer, same units as the
  corresponding `leafbio` entries).

- spectral:

  list. Spectral configuration, including wavelength vectors `wlP`,
  `wlE`, `wlF`, `wlS` and region boundaries `reg1`/`reg2`/`reg3`/`IwlT`.

- leafbio:

  list. Baseline leaf biochemical/structural properties; the fields
  `Cab`, `EWT`, `Car`, `LMA`, `Cs`, `N` are overwritten per layer from
  `mly`.

- soil:

  list. Soil properties; only `rs_thermal` (thermal-region soil
  reflectance) is used, and only for the optional diagnostic plots.

- optipar:

  list. Leaf-level optical parameters used by the underlying Fluspect
  model.

- nl:

  integer. Number of canopy layers over which the per-profile-layer leaf
  optics are replicated.

- step:

  numeric, optional. Wavelength step (nm) used to compute the `Mb`/`Mf`
  excitation-emission matrices; if missing, a default step of 5 nm is
  used (matrix of 53 x 71), a step of 1 nm gives a 211 x 351 matrix.
  Must be lower than 8 nm.

- get.plots:

  logical. If `TRUE`, produces diagnostic plots of leaf reflectance,
  pigment contribution factors, PSI/PSII quantum yield fractions and the
  `Mb`/`Mf` matrices. Default `TRUE`.

## Value

A list `leafopt` with: `refl`/`tran` (leaf reflectance/transmittance per
canopy layer x wavelength, `[nl x length(wlP)]`), `kChlrel`/`kCarrel`
(relative contribution of chlorophyll/carotenoids to absorption, same
dimensions), `Mb`/`Mf` (backward/forward fluorescence
excitation-emission matrices, replicated per canopy layer),
`phiI`/`phiII` (PSI/PSII relative quantum yield spectra at the
fluorescence wavelengths).

## Author

    Christiaan van der Tol  (Original version in Matlab)

Carlos Camino (Ported version into R)

## Examples

``` r
if (FALSE) { # \dontrun{
leafopt <- get.fluspect_mSCOPE(mly, spectral, leafbio, soil, optipar, nl,
                                step = 5, get.plots = FALSE)
} # }
```
