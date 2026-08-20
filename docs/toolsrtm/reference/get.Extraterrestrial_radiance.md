# get.Extraterrestrial_radiance

get.Extraterrestrial_radiance

## Usage

``` r
get.Extraterrestrial_radiance(Ea0, DOY, tts)
```

## Arguments

- Ea0:

  Solar constant spectral irradiance

- DOY:

  day of year integer.

- tts:

  solar zenith angle

## Value

solar extraterrestrial spectrum

## Examples

``` r
# Example usage:
# Assuming Ea0, DOY, and tts are your input values
Ea <- get.Extraterrestrial_radiance(Ea0 = 1361, DOY = 180, tts = 45)
print(Ea)
#> [1] 296.1273
# 296.1273
```
