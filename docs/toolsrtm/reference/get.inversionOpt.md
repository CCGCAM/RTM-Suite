# LUT Inversion Using a Radiative Transfer Model (RTM)

This function performs the inversion of a Radiative Transfer Model (RTM)
based on observed sensor reflectance values. It compares simulated
reflectance values from the RTM with observed values and selects the
best matches using different merit functions.

## Usage

``` r
get.inversionOpt(
  rfl.sensor = NULL,
  rfl.rtm = NULL,
  LUT = NULL,
  wave = NULL,
  method = "merit-RMSE",
  nOpt = NULL,
  custom_stat = NULL
)
```

## Arguments

- rfl.sensor:

  A matrix with reflectance values of the observed sensor (e.g., rows
  representing different observations, columns representing
  wavelengths).

- rfl.rtm:

  A matrix with reflectance values simulated by the Radiative Transfer
  (RT) model. The structure of this matrix should be similar to
  `rfl.sensor`.

- LUT:

  A LUT (Look-Up Table) containing the distribution of biophysical
  parameters used as input in the RT model (e.g., leaf chlorophyll
  content, water content).

- wave:

  A vector containing the wavelengths corresponding to the columns of
  `rfl.sensor` and `rfl.rtm`.

- method:

  The merit function used to evaluate the inversion. Options include:

  - `'merit-RMSE'`: Root Mean Square Error. The default when no method
    is provided

  - `'merit-NRMSE'`: Normalized RMSE (scaled by the range of observed
    data).

  - `'merit-MAE'`: Mean Absolute Error.

  - `'merit-NMB'`: Normalized Mean Bias.

  - `'merit-FGE'`: Fractional Gross Error.

  - `'merit-DWT'`: RMSE computed on discrete wavelet transform (Haar)
    coefficients instead of raw reflectance.

  - `'merit-1stD'`: RMSE computed on the first derivative of the spectra
    (finite differences along `wave`) instead of raw reflectance.

  - `'merit-custom.metric'`: A custom metric defined by users.

- nOpt:

  The number of optimal solutions (i.e., the best-matching simulated
  spectra) to select based on the chosen merit function.

- custom_stat:

  An optional custom statistic function. If provided, this will override
  the default merit function. The custom function should take two
  arguments: the simulated and observed values and return a single
  numeric value (the error or difference metric).

## Value

A list with two elements:

- `rfl.b`: A matrix of the best-matching reflectance values selected
  from `rfl.rtm`.

- `LUT.best`: A data frame containing the corresponding biophysical
  parameters from the LUT for the best solutions.

## Examples

``` r
# Simulated example usage:
sensor_data <- matrix(runif(100), nrow = 10, ncol = 10) # Simulated sensor reflectance
rtm_data <- matrix(runif(100), nrow = 10, ncol = 10)    # Simulated RTM reflectance
lut_table <- data.frame(N = runif(10), Cab = runif(10), Cw = runif(10)) # Simulated LUT
wavelengths <- seq(400, 700, length.out = 10)  # Simulated wavelengths
result <- get.inversionOpt(sensor_data, rtm_data, lut_table, wavelengths,
                            method = 'merit-RMSE', nOpt = 5)
#> Merit function using merit-RMSE is processing
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======                                                               |  10%  |                                                                              |==============                                                        |  20%  |                                                                              |=====================                                                 |  30%  |                                                                              |============================                                          |  40%  |                                                                              |===================================                                   |  50%  |                                                                              |==========================================                            |  60%  |                                                                              |=================================================                     |  70%  |                                                                              |========================================================              |  80%  |                                                                              |===============================================================       |  90%  |                                                                              |======================================================================| 100%
print(result)
#> [[1]]
#>           R.400 R.433.333333333333 R.466.666666666667     R.500
#>  [1,] 0.5777439          0.4379375          0.5403002 0.5479565
#>  [2,] 0.3388092          0.5918006          0.5705170 0.4771207
#>  [3,] 0.6814820          0.4227113          0.5941492 0.4818942
#>  [4,] 0.5166098          0.6318400          0.6464851 0.4454986
#>  [5,] 0.3420301          0.7026746          0.6239810 0.3885904
#>  [6,] 0.5208623          0.6535061          0.5967224 0.3470724
#>  [7,] 0.7070879          0.5225232          0.4870611 0.5008718
#>  [8,] 0.2358851          0.4409694          0.6524779 0.5063533
#>  [9,] 0.5250348          0.4608177          0.4608557 0.6309201
#> [10,] 0.3420301          0.7026746          0.6239810 0.3885904
#>       R.533.333333333333 R.566.666666666667     R.600 R.633.333333333333
#>  [1,]          0.3828930          0.5275253 0.3758500          0.3341890
#>  [2,]          0.5427068          0.2695723 0.4814252          0.5032734
#>  [3,]          0.4472590          0.5150185 0.3506139          0.3783419
#>  [4,]          0.5663552          0.4989844 0.6158093          0.3336079
#>  [5,]          0.5763804          0.4488285 0.4855131          0.4104537
#>  [6,]          0.5356454          0.3458877 0.5227670          0.4754391
#>  [7,]          0.3796069          0.4438090 0.5782080          0.3613078
#>  [8,]          0.6229248          0.4983867 0.6880725          0.4249039
#>  [9,]          0.3866682          0.3674936 0.5368662          0.3891422
#> [10,]          0.5763804          0.4488285 0.4855131          0.4104537
#>       R.666.666666666667     R.700
#>  [1,]          0.5085053 0.2085346
#>  [2,]          0.3933704 0.3826236
#>  [3,]          0.5845443 0.3776429
#>  [4,]          0.7612825 0.3485879
#>  [5,]          0.7334558 0.2451866
#>  [6,]          0.5287351 0.5569679
#>  [7,]          0.4319006 0.4994925
#>  [8,]          0.5010785 0.3498130
#>  [9,]          0.2965359 0.3251481
#> [10,]          0.7334558 0.2451866
#> 
#> [[2]]
#>    ID_lut merit-RMSE     N   Cab    Cw
#> 1     5.4      0.368 0.647 0.548 0.397
#> 2     4.6      0.359 0.629 0.515 0.513
#> 3     5.0      0.367 0.550 0.399 0.403
#> 4     6.2      0.392 0.653 0.353 0.386
#> 5     7.2      0.348 0.703 0.381 0.506
#> 6     5.0      0.334 0.492 0.345 0.360
#> 7     4.2      0.384 0.482 0.500 0.302
#> 8     5.6      0.295 0.653 0.516 0.630
#> 9     3.8      0.387 0.619 0.669 0.455
#> 10    7.2      0.342 0.703 0.381 0.506
#> 
```
