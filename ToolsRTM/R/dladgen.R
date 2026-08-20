#' Computes the leaf angle distribution function value (freq)
#'
#' Using the original bimodal distribution function initially proposed in SAIL
#'  References
#'  ----------
#'  (Verhoef1998) Verhoef, Wout. Theory of radiative transfer models applied
#'  in optical remote sensing of vegetation canopies.
#'  Nationaal Lucht en Ruimtevaartlaboratorium, 1998.
#'  http://library.wur.nl/WebQuery/clc/945481.
#' @param a controls the average leaf slope
#' @param b controls the distribution's bimodality
#' LIDF type 		  a 		b
#' Planophile 	  1		  0
#' Erectophile    -1	 	0
#' Plagiophile 	  0		  -1
#' Extremophile 	0		  1
#' Spherical 	    -0.35 -0.15
#' Uniform        0     0
#' requirement: |LIDFa| + |LIDFb| < 1
#'
#' @return foliar_distrib list. lidf and litab
#' @export
dladgen <- function(a,b){

   # a=1
  #  b=2

litab <- c(5,15,25,35,45,55,65,75,81,83,85,87,89)
freq=c() ##added
for (i1 in c(1:8)){
    t   <- i1 * 10
    freq[i1] <- ToolsRTM::dcum(a,b,t)
}

for (i2 in c(9:12)){
    t   <- 80 + (i2 - 8) * 2
    freq[i2] <- ToolsRTM::dcum(a,b,t)
}

freq[13] <- 1
for (i   in seq(13,2,-1)){
    freq[i] <- freq[i] - freq[i-1]
}


LeafDistribution <- list("lidf" = freq,"litab" =litab) ##added
return(LeafDistribution) ##added

}

