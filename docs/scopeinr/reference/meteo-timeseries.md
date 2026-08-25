# Example meteorological/radiation time series for diurnal SCOPE simulations

Example time series (one value per time step) used by
[`getLUT_time`](getLUT_time.md) to drive multi-timestep SCOPE runs:
incoming shortwave radiation (`Rin_`), incoming longwave radiation
(`Rli_`), sun/sky irradiance spectra (`Esun_`/`Esky_`), air temperature
(`Ta_`), vapour pressure (`ea_`), air pressure (`p_`), wind speed
(`u_`), time-of-day (`t_`) and year (`year_`).

## Usage

``` r
Rin_

Rli_

Esun_

Esky_

Ta_

ea_

p_

u_

t_

year_
```

## Format

A numeric vector or data frame, one row/element per time step.

An object of class `data.frame` with 214 rows and 1 columns.

An object of class `data.frame` with 2162 rows and 1 columns.

An object of class `data.frame` with 2162 rows and 1 columns.

An object of class `data.frame` with 214 rows and 1 columns.

An object of class `data.frame` with 214 rows and 1 columns.

An object of class `data.frame` with 214 rows and 1 columns.

An object of class `data.frame` with 214 rows and 1 columns.

An object of class `data.frame` with 214 rows and 1 columns.

An object of class `data.frame` with 216 rows and 1 columns.
