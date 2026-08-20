#' get.spectra.SCOPE
#' \code{get.spectra.SCOPE} Calculates the spectra of hemisperical and directional observed #' a function to get the spectral characteristics for SCOPE model
#'
#' @param getSpectral  is T get the spectral object needed for SCOPE
#'
#' @return a list with al spectral ranges and Intervals
#' @export
#'
#' @examples
#' spectral <- get.spectra.SCOPE(getSpectral = TRUE)
get.spectra.SCOPE<-function(getSpectral=T){

  if (getSpectral == T){

    #Define spectral regions for SCOPE v_1.40
    # All spectral regions are defined here as row vectors
    # WV Jan. 2013
    spectral<-list()
    # 3 spectral regions for SCOPE

    ### region 1
    step1=1 ## 1nm
    reg1<-seq(400,2400,step1)
    spectral$reg1<-reg1
    ### region 2
    step2=100 ## 100nm
    reg2<- seq(2500,15000,step2)
    spectral$reg2<-reg2

    ### region 3
    step3=1000 ## 100nm
    reg3<-seq(16000,50000,step3)
    spectral$reg3<-reg3

    spectral$wlS<-c(reg1, reg2, reg3)

    # Other spectral (sub)regions
    spectral$wlP   = reg1                            # PROSPECT data range
    spectral$wlE   = seq(400,750,1)                       # excitation in E-F matrix
    spectral$wlF   = seq(640,850,1)                       # chlorophyll fluorescence in E-F matrix
    spectral$wlO   = reg1                           # optical part
    spectral$wlT   = c(reg2, reg3)                     # thermal part
    spectral$wlZ   = c(500,600,1)                       # xanthophyll region
    wlS            = spectral$wlS;
    spectral$wlPAR = wlS[wlS>=400 & wlS<=700];  # PAR range
    spectral$IwlF  = c(640:850)-399
    spectral$IwlP  = 1:length(reg1);

    # length of spectral$IwlT is the number of elements in the concatenation of reg2 and reg3.
    #spectral$IwlT  = c(length(reg1)+1:length(reg1)+length(reg2)+length(reg3))
    spectral$IwlT <- (length(reg1) + 1):(length(reg1) + length(reg2) + length(reg3))


    wave=c(seq(from=400, to=2400,by=1),seq(from=2500, to=15000,by=100),seq(from=16000, to=50000,by=1000))
    spectral$wlIrrad =wave
    # Data used by aggreg routine to read in MODTRAN data

    spectral$SCOPEspec.nreg = 3
    spectral$SCOPEspec.start = c(400,  2500,  16000)
    spectral$SCOPEspec.end   = c(2400, 15000,  50000)
    spectral$SCOPEspec.res   = c(1,100,   1000)

    return(spectral)
  }

}
