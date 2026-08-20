# Introduction to InversionOpt

## Working with InversionOpt

\*\* InversionOpt \*\*

Is a LUT-based method for estimating the best parameters from PROSAIL
model using three methods:

``` r

library(ToolsRTM)

##method is the method (opt ='merit-RMSE','merit-DWT',merit-1stD')

#inv.RMSE<-InversionOpt_nOpt(rfl.sensor=rfl.sensor, #observado
 #                           rfl.prosail=rfl.prosail, #500
  #                          LUT=LUT,
   #                         wave=wave.vnir.swir, 
    #                        n=n_cases, #500
     #                       method='merit-RMSE', 
        #                    nOpt=100)
```
