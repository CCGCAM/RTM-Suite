1. Source of errors for SCOPEinR:

 - 1) SCOPEinR::get.fluspect_mSCOPE  --> fixed error for Mb and Mf, should return an array of Mf[53,71, nl], each nl is
                                         a repetition of 53x71. 53x71 is moduled by step=5 (interpolation each 5 wave)
                                         (SOLVED, Claude Code session 2026-08: the array only kept the LAST mSCOPE
                                         profile layer's Mb/Mf for the whole canopy instead of each layer's own,
                                         replicated across its own sublayers -- fixed to match the rho/tau pattern.)

      Here ..... I found and error in Mf and Mb, change values in Phi to  optipar to get similar result as SCOPE (dafault)
      optipar2017.ProspectD$phi, the errors are in mfluspect_mSCOPE
      (Phi/optipar version question NOT resolved -- still an open physics/parameterization call for Carlos, see below.)


 - 2) SCOPEinR::meanleaf  --> getRTMo and get.ebal should be verify, the three possiblilities
 - 3) SCOPEinR::get.aggregator.ebal() --> get.ebal, check the final results, dependencies with meanleaf
 - 4) SCOPEinR::get.biochemical() --> check the implementation of brent fucntion, dependencies
      in get.ebal
 - 5) SCOPEinR::get.biochemical() -- > check eta for data.bcu, bca (solved)

 - 6) SCOPEinR::get.RTMf() -- > check values in sfEs and sbEs, both equation use same inputs (solved)
 - 7) SCOPEinR::get.RTMt.planck() --> check the Hcsu value ii is a matrix (section 1.2), possible source of errors

 - 8) SCOPEinR::get.ebal() --> Energy balance error soil is huge. Possibilities is in resistance using get.resistance
                           --> Fixed maxEBer_soil = 200 to avoid printing in SCOPE model, please check what happen,                                      possibiliites RTMo and ebalRemoved options
                             To avoid to sho message I change  if (counter >= maxit)  by if ((counter <= maxit) ), Here
                             I foind errors in the WHile con, always give a counter higer tham maxit
                             - Alternatives in options.file
                             (SOLVED, Claude Code session 2026-08: root cause found. The soil temperature update
                             `Ts <- Ts[1] + Wc*EBers/(...)` only ever used Ts[1] (shaded) in the numerator for BOTH
                             the shaded and sunlit components, while the denominator correctly used the full Ts
                             vector -- so the sunlit-soil temperature never updated from its own previous value.
                             Fixed to `Ts <- Ts + Wc*EBers/(...)`, matching how Tch/Tcu are updated. Default
                             get.SCOPE() run now converges in counter=7 iterations, maxEBers=0.108 W/m^2 (was
                             counter=101/never-converged, maxEBers=152.13 W/m^2). The counter>=maxit print
                             condition was also flipped back to correct in a separate fix.)

#### Errors in get.outs values
 - 1) SigmaF -->   scape probabilities  check , the values is so much difference !!!!
 - 2) df.scalers -->   Values so low for all values, F_1stPeak wl_1stPeak return NA values for 688 F
                       Main codes in SCOPE.m and RTMf
 - 3) df.LoF --> Values of fluoresecen so low, check RMTf values for LoF_, Femliave_, EoutF_, EoutFrc_


   possibildades estan en Mf and MB, las primeras dimesiones no estan bien,


 - check fluoresecence associate to phI, and PhiII (Fluospect, selected optipar )
