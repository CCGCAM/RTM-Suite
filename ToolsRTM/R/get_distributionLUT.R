#' Build a LUT with a per-trait distribution choice (Uniform or Gaussian) and an optional Car~Cab correlation
#'
#' Alternative to \code{\link{getLUT}} for when different traits need
#' different sampling distributions in the same LUT (e.g. LAI sampled
#' Uniform while Cab is sampled Gaussian), rather than one distribution for
#' every trait. Gaussian traits are drawn by truncated-normal rejection
#' sampling (\code{\link{gauss_byMin_Max}}), so they still respect
#' \code{minval}/\code{maxval} bounds.
#'
#' @param minval one-row data.frame/list, column names = trait names, minimum value per trait.
#' @param maxval one-row data.frame/list, column names = trait names, maximum value per trait.
#' @param nSamples numeric. Number of LUT rows to generate.
#' @param TypeDistrib named list, one entry per trait in \code{minval}, each either
#'   \code{"Uniform"} or \code{"Gaussian"}.
#' @param Mean_gauss one-row data.frame/list, mean per trait -- only read for traits where
#'   \code{TypeDistrib} is \code{"Gaussian"}.
#' @param Std_gauss one-row data.frame/list, standard deviation per trait -- only read for
#'   traits where \code{TypeDistrib} is \code{"Gaussian"}.
#' @param DepCab logical. If \code{TRUE} and \code{"Car"} is one of the traits, Car is not
#'   drawn independently -- it's redrawn as \code{\link{correlatedValue}(Cab/4, r = 0.8)},
#'   the empirical Cab-Car co-variation seen in leaf pigment data.
#' @param setseed integer. Random seed.
#'
#' @return LUT as a data.frame, one column per trait in \code{minval}.
#' @importFrom stats runif rnorm sd
#'
#' @export
#'
#' @examples
#' minv <- data.frame(Cab = 10, Car = 2, LAI = 0.5)
#' maxv <- data.frame(Cab = 80, Car = 20, LAI = 7)
#' distrib <- list(Cab = "Gaussian", Car = "Uniform", LAI = "Uniform")
#' meang <- data.frame(Cab = 40, Car = NA, LAI = NA)
#' stdg  <- data.frame(Cab = 15, Car = NA, LAI = NA)
#' LUT <- get_distributionLUT(minval = minv, maxval = maxv, nSamples = 100,
#'                             TypeDistrib = distrib, Mean_gauss = meang, Std_gauss = stdg,
#'                             DepCab = TRUE, setseed = 1)
get_distributionLUT<-function(minval=NULL,maxval=NULL,nSamples=NULL,TypeDistrib=NULL,Mean_gauss=NULL, Std_gauss=NULL, DepCab=NULL, setseed=NULL){
  # define InputPROSAIL # 3 random parameters
   inputLUT<-list()
   
   n_casesNorm=nSamples*2
   
   value_seed=setseed
   set.seed(value_seed)
 
   for (i in 1:length(minval)){
      
      trait <- names(minval)[i]
     
      # if uniform distribution
      if(TypeDistrib[[trait]] == 'Uniform') {
   
        inputLUT[[trait]] <- stats::runif(nSamples,min = minval[1,trait],max=maxval[1,trait])
         if (names(inputLUT)[i] == 'Car' & DepCab == T){
        
         inputLUT[[trait]] <- ToolsRTM::correlatedValue(x=inputLUT[['Cab']]/4, r=.8) 
         }
        
       
      }
      # if Gaussian distribution
      else  {
   
      inputLUT[[trait]] <- ToolsRTM::gauss_byMin_Max(n=nSamples, m=Mean_gauss[1,trait], s=Std_gauss[1,trait], lwr=minval[1,trait], upr=maxval[1,trait], nnorm=n_casesNorm)
      if (names(inputLUT)[i] == 'Car' & DepCab == T){
         set.seed(value_seed)
         inputLUT[[trait]] <- ToolsRTM::correlatedValue(x=inputLUT[['Cab']]/4, r=.8) 
      }
      
      }
      
    
   }
   LUT.dataframe <-data.frame(do.call(cbind,inputLUT))
   set.seed(Sys.time())
  return(LUT.dataframe)
}