# Compute vegetation indices from reflectance spectra (v2, expanded index set)

Interpolates reflectance to 1nm resolution and computes a large set of
published vegetation indices (VNIR and/or SWIR domain), one row per
input spectrum. This is an expanded version of
[`getIndices()`](getIndices.md) with additional SWIR-domain and red-edge
indices.

## Usage

``` r
get.indices.v2(data, pattern.rfl = "R.", factor = NULL, spectral.domain = NULL)
```

## Arguments

- data:

  data.frame or matrix of reflectance spectra, one row per sample, with
  reflectance columns named by wavelength (e.g. "R.400", "R.405", ...).

- pattern.rfl:

  character. Prefix identifying reflectance columns in `data`. Default
  "R.".

- factor:

  numeric. Multiplier applied to reflectance values before computing
  indices (e.g. to rescale 0-10000 integer reflectance to 0-1). Default
  1 (no rescaling).

- spectral.domain:

  character. One of "VNIR" (400-850nm), "SWIR" (800-2550nm), or
  "VNIR-SWIR" (400-2550nm). Default "VNIR".

## Value

A data.frame of computed indices, one row per input spectrum, with
columns that are entirely NA removed.

## Examples

``` r
# Synthetic reflectance spectrum, VNIR domain
wl <- seq(400, 850, by = 5)
rfl <- matrix(runif(length(wl), 0.05, 0.5), nrow = 1)
colnames(rfl) <- paste0("R.", wl)
idx <- get.indices.v2(as.data.frame(rfl), pattern.rfl = "R.", spectral.domain = "VNIR")
#> [1] "Estimating indices using VNIR domain..."
#>   |                                                                              |                                                                      |   0%  |                                                                              |======================================================================| 100%
#>            NDVI       RDVI        SR       MSR      OSAVI      MSAVI     MTVI1
#> [1,] -0.2829068 -0.2379507 0.5589596 -0.252364 -0.2676401 -0.2303549 -0.546953
#>           MTVI2      MCARI     MCARI1     MCARI2       EVI      LIC1      VOG
#> [1,] -0.3594276 -0.1041945 -0.6836913 -0.3594276 -1.531349 0.4040859 1.093503
#>           VOG2      VOG3       GM1       GM2      TCARI      T.O        CI
#> [1,] 0.1770907 0.1996262 0.3015598 0.2946145 -0.4798857 1.793026 0.2110397
#>            TVI      SRPI       NPQI      NPCI     CTR1      CAR   DCabxc
#> [1,] -28.65244 0.9707103 -0.1179667 0.0148625 1.533203 1.002605 1.606571
#>       DNCabxc      SIPI     CRI550     CRI700  CRI550m CRI700m  RCRI550
#> [1,] 14.35092 -0.198695 -0.2520837 -0.1720908 8.768547 8.84854 1.553493
#>       RCRI700        PSRI     LIC3      CIre CIrededge CIgreen Chlred.edge
#> [1,] 1.591901 -0.08961966 1.451741 0.7207813 0.5972911      NA   0.6092661
#>           CVI      IRECI      REP      RVI    RedEg1     RedEg2        PRI
#> [1,] 1.030785 0.02340351 739.2644 1.030785 0.6453334 -0.2155591 -0.6092494
#>          PRI515      PRIM1      PRIM2     PRIM3      PRIM4     PRIn    PRI_CI
#> [1,] -0.6084309 -0.1408465 0.07562922 0.2190063 -0.8027216 3.942487 0.3880476
#>              B         G         R    BGI1      BGI2      BF1      BF2      BF3
#> [1,] 0.2855239 0.6344809 0.6494384 1.69012 0.2274706 6.468157 2.665213 4.656677
#>           BF4      BF5     BRI1      BRI2       RGI      RARS     LIC2
#> [1,] 1.313835 7.430059 4.503333 0.6060965 0.3753043 0.7885633 3.427625
#>              HI     CUR    PSSRa     PSSRb     PSSRc       PSNDc CR.red.nir.1
#> [1,] -0.3906864 0.89648 2.356188 0.6218097 0.8976747 -0.05392143    0.8964065
#>      CR.red.nir.2 CR.red.nir.3 CR.red.nir.4 CR.red.nir.5 CR.red.nir
#> [1,]    0.8899571    0.7264426     2.426856       1.2071   1.134038
#>      CR.red.nir.7
#> [1,]    0.9581846
```
