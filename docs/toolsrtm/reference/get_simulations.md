# Simulate spectra with RT models

Simulate spectra with RT models

## Usage

``` r
get_simulations(inputLUT = NULL, psoil = 0.5, rtm.model = "PROSAIL")
```

## Arguments

- inputLUT:

  is a LUT with all the inputs needed for simulating a leaf or canopy
  model

- psoil:

  soil factor, by default is set to 0.5

- rtm.model:

  model for generating the simulations ('PROSPECT-PRO','PROSAIL' or
  'INFORM')

## Value

a LUT with all the simulations and the plot
