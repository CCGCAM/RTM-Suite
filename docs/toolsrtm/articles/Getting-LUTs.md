# Getting-LUTs

Hellloss

``` r

library(ToolsRTM)
```

define min and max values for all parameters defined in TypeDistrib

``` r

minval <- data.frame('N' = 1.5,'Cab'=5,'Car'=0,'Anth' = 0,'Cbrown'= 0.0,
                     'EWT' = 0.001,'Prot' =  0.00001, 'CBC' = 0.00001,
                     'LIDFa' = 0, 'LAI' = 0.5,
                     ### input for INFORM
                     LAIu = 0, sd = 200,cd = 0.2 , h=5,
                     ### input for fourSAIL-2
                     fraction_brown = 0, diss = 0.1,Cv = 0.3, Zeta=0)
```

define min and max values for all parameters defined in TypeDistrib

``` r

maxval <- data.frame('N' = 3,'Cab'=70,'Car'=25,'Anth' = 7,'Cbrown'= 0.2,
                     'EWT' = 0.035,'Prot' =  0.03, 'CBC' = 0.03,
                     'LIDFa' = 70, 'LAI' = 7,
                      ### input for INFORM
                      LAIu = 0.8, sd = 1000,cd = 7, h=20,
                      ### input for fourSAIL-2
                     fraction_brown = 1, diss = 1,Cv = 1, Zeta=0.2)
```

## define the type of distribution

``` r

TypeDistrib<-data.frame('N' = 'Gaussian','Cab'='Gaussian','Car'='Gaussian',
                        'Anth' = 'Uniform','Cbrown'= 'Uniform',
                        'EWT' = 'Uniform',
                        'Prot' =  'Uniform', 'CBC' = 'Uniform',
                        'LIDFa' = 'Uniform', 'LAI' = 'Gaussian',
                        ### input for INFORM
                        LAIu = 'Uniform', sd = 'Uniform',cd = 'Uniform' , h='Uniform',
                        ### input for fourSAIL-2
                        fraction_brown = 'Uniform', diss = 'Uniform',Cv = 'Uniform', Zeta='Uniform')
```

define mean and STD for gaussian distributions

``` r

Mean_gauss <- data.frame('N'=2.2,'Cab'=45,'Car'=8,'LAI' = 2.25)
std_gauss <- Mean_gauss/2.0
```

``` r

nSamples=500
j=2
data.LUT<-get_distributionLUT(minval=minval,maxval=maxval,
                              nSamples=nSamples,TypeDistrib=TypeDistrib,
                              Mean_gauss=Mean_gauss, Std_gauss=std_gauss,DepCab = T,setseed = j*123)
```
