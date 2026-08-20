# get.Fluorescence.model `get.Fluorescence.model` Fluorescence model

get.Fluorescence.model `get.Fluorescence.model` Fluorescence model

## Usage

``` r
get.Fluorescence.model(ps, x, Kp, Kf, Kd, Knparams)
```

## Arguments

- ps:

  numeric. Photochemical yield (fraction of absorbed light used in
  photochemistry).

- x:

  numeric. Normalized (0-1) intermediate variable driving
  non-photochemical quenching, passed in to avoid recomputing it
  internally.

- Kp:

  numeric. Rate constant for photochemistry.

- Kf:

  numeric. Rate constant for fluorescence.

- Kd:

  numeric. Rate constant for (light-independent) thermal dissipation.

- Knparams:

  numeric vector of length 3: Kno (max NPQ rate constant), alpha, beta -
  shape parameters of the non-photochemical quenching response.

## Value

A list with eta, qE, qQ, fs, fo, fm, fo0, fm0, and Kn - fluorescence
yields and quenching parameters.

## Examples

``` r
get.Fluorescence.model(ps = 0.3, x = 0.5, Kp = 4, Kf = 0.05, Kd = 0.95,
                        Knparams = c(2.48, 2.83, 0.114))
#> $eta
#> [1] 1.385679
#> 
#> $qE
#> [1] 0.696661
#> 
#> $qQ
#> [1] 0.4894379
#> 
#> $fs
#> [1] 0.01385679
#> 
#> $fo
#> [1] 0.007661851
#> 
#> $fm
#> [1] 0.01979541
#> 
#> $fo0
#> [1] 0.01
#> 
#> $fm0
#> [1] 0.05
#> 
#> $Kn
#> [1] 1.525838
#> 
```
