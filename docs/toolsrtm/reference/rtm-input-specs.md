# RTM parameter ranges and prior distributions (leaf/canopy models, shared)

`inputsFlUSPECT`, `inputsINFORM`, `inputsLiberty`, `inputsPROSAIL`,
`inputsSPART` and `inputsRTMs` each give the default parameter bounds
and sampling distributions used by [`getLUT()`](getLUT.md) and related
look-up-table generators for the corresponding model (Fluspect, INFORM,
Liberty, PROSAIL, SPART, and a combined table across models
respectively).

## Usage

``` r
inputsFlUSPECT

inputsINFORM

inputsLiberty

inputsPROSAIL

inputsSPART

inputsRTMs
```

## Format

A data frame, one row per model parameter, with columns `variable`
(parameter name), `lower`/`upper` (sampling bounds), `units`,
`Distribution` (sampling distribution family), `Mean_D`/`Std_D`
(distribution mean/SD where applicable), `Dependencies` (other
parameters this one covaries with, if any), `use.default` (whether the
package default is used) and `default` (default value). `inputsRTMs`
additionally has a `model` column identifying which RTM each row belongs
to.

An object of class `data.frame` with 24 rows and 10 columns.

An object of class `data.frame` with 25 rows and 10 columns.

An object of class `data.frame` with 18 rows and 10 columns.

An object of class `data.frame` with 30 rows and 10 columns.

An object of class `data.frame` with 40 rows and 11 columns.
