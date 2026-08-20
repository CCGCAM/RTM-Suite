
#' Simulate reflectance spectra while sweeping one trait, all others fixed
#'
#' Builds a LUT where every trait is held at its default except \code{trait},
#' which is swept from \code{nmin} to \code{nmax} in steps of \code{Interval},
#' runs the canopy simulation, and plots the resulting spectra.
#'
#' @param trait character. Name of the trait to sweep (e.g. \code{"Cab"}, \code{"LAI"}). Default \code{"Cab"}.
#' @param nmin,nmax numeric. Minimum/maximum value of \code{trait}'s sweep range.
#' @param Interval numeric. Step size between successive values of \code{trait} in the sweep.
#' @param psoil numeric (0-1). Soil brightness factor (0 = wet soil, 1 = dry soil), mixed
#'   from \code{ToolsRTM::dataSpec_PDB}'s dry/wet reference spectra.
#' @param model character. Canopy simulation to run: \code{"PROSAIL"}, \code{"INFORM"}, or \code{"PROSPECT"}.
#' @param method character. \code{"ggplot"} (default) or \code{"classical"} (base-R plotting) for the spectral plot.
#' @export
#'
#' @return If \code{method = "ggplot"}: a list with \code{LUT} (the swept-trait LUT used) and
#'   \code{Plot} (the ggplot object). If \code{method = "classical"}: just the LUT data.frame
#'   (the plot is drawn as a side effect via base graphics, not returned).
#'
#' @examples
#' \dontrun{
#' out <- getSim_fromLUT(trait = "Cab", nmin = 10, nmax = 80, Interval = 10, model = "PROSAIL")
#' out$Plot
#' }
getSim_fromLUT<-function(trait='Cab',nmin=0, nmax=100, Interval=10, psoil=0.5,model='PROSAIL', method='ggplot'){
  
  
  require("foreach")
  
  data <- ToolsRTM::dataSpec_PDB
  Rsoil1  <- data[,11]  # rsoil1 = dry soil
  Rsoil2 <- data[,12]  # rsoil2 = wet soil 
  
  if ( is.null(psoil)){
    message(' soil factor of 0.5 will be used as default')
    psoil	 <-  0.5    # soil factor (psoil=0: wet soil / psoil=1: dry soil)
    rsoil  <- psoil*Rsoil1+(1-psoil)*Rsoil2
  } else {
    rsoil  <- psoil*Rsoil1+(1-psoil)*Rsoil2
  }
  
  set.seed(12345)
  
  if ( is.null(model)){
    message(' PROSAIL model will be used as input')
  }
  
  #Inputs by default for Leaf model (PROSPECT)
  N=2.5;  Cab =40; Car=15; Anth=0; Cbrown=0.5
  EWT=0.09; LMA=0.012; alpha=40
  Prot=0;CBC=0
  #Inputs by default for canopy model (fourSAIL) 
  ## parameters for fourSAIL
  LIDFa=30; LIDFb=0; TypeLidf=2; lai=3
  hotspot=0; tts=5; tto=0; psi=0
  
  ## Input for INform model
  LAIu=0.5
  sd=500; cd=4.5; h=20; psi=30
  skyl=0.1
  
  
  if (is.null(Interval)){
    message('number of classes is null, by default 10 classes will be greated')
    
    Interval = 10
    nLUT = 1 + Interval
  }
  if ( is.null(trait)){
    message(' trait is empty, chorophyll content will be used as default')
    input<- seq(nmin,nmax, by=Interval)##stats::runif(nLUT,min = nmin,max=nmax)
    Cab=input
  } else if (trait == 'Cab'){
    input<- seq(nmin,nmax, by=Interval)##stats::runif(nLUT,min = nmin,max=nmax)
    Cab=input
  } else if (trait == 'N'){
      input<- seq(nmin,nmax, by=Interval)
      N=input
  } else if (trait == 'Car'){
      input<- seq(nmin,nmax, by=Interval)#stats::runif(nLUT,min = nmin,max=nmax)
      Car=input
  } else if (trait == 'Anth'){
    input<- seq(nmin,nmax, by=Interval)
    Anth=input
  } else if (trait == 'Cbrown'){
    input<- seq(nmin,nmax, by=Interval)
    Cbrown=input
  } else if (trait == 'EWT'){
    input<-seq(nmin,nmax, by=Interval)
    EWT=input
  } else if (trait == 'LMA'){
    input<- seq(nmin,nmax, by=Interval)
    LMA=input
  } else if (trait == 'alpha'){
    input<- seq(nmin,nmax, by=Interval)
    alpha=input
  } else if (trait == 'Prot'){
    LMA=0
    CBC= 0.02
    input<- seq(nmin,nmax, by=Interval)
    Prot=input
  } else if (trait == 'CBC'){
    LMA=0
    Prot = 0.02
    input<- seq(nmin,nmax, by=Interval)
    CBC=input
  } else if (trait == 'LAI'){
    input<- seq(nmin,nmax, by=Interval)
    lai=input
  } else if (trait == 'LIDFa'){
    input<- seq(nmin,nmax, by=Interval)
    LIDFa=input
  } else if (trait == 'hotspot'){
    input<- seq(nmin,nmax, by=Interval)
    hotspot=input
  } else if (trait == 'tts'){
    input<- seq(nmin,nmax, by=Interval)
    tts=input
  } else if (trait == 'tto'){
    input<- seq(nmin,nmax, by=Interval)
    tto=input
  } else if (trait == 'psi'){
    input<- seq(nmin,nmax, by=Interval)
    psi=input
  } else if (trait == 'LIDFb'){
    message('Type LIDF is set to 2; and LIDFb = 0')
    LIDFb = 0
  } else if (trait == 'LAIu'){
    input<- seq(nmin,nmax, by=Interval)
    LAIU=input
  } else if (trait == 'sd'){
    input<- seq(nmin,nmax, by=Interval)
    sd=input
  } else if (trait == 'cd'){
    input<- seq(nmin,nmax, by=Interval)
    cd=input
  } else if (trait == 'h'){
    input<- seq(nmin,nmax, by=Interval)
    h=input
  } 
  
  
  if ( is.null(model) | model == 'PROSAIL'){
    message('By default, PROPECT + fourSAIL model is processing ...')
    LUT<- data.frame(N,Cab,Car,Anth,Cbrown,EWT,LMA,alpha,Prot,CBC,
                     LIDFa, LIDFb=0, TypeLidf=2, LAI=lai,
                     hspot=hotspot, tts, tto, psi)
    nLUT = dim(LUT)[1]
    ## choose number of processors/cores
    no_cores <- parallel::detectCores() - 2 
    cl <- parallel::makeCluster(no_cores)
    doParallel::registerDoParallel(cl)
    
    start_time <- Sys.time()
    sim.rfl<-list()
    ### simulations
      sims<-foreach(i=1:nLUT) %dopar% {
        data.foursail_pro<-ToolsRTM::foursail(inputLUT=LUT[i,],rsoil=rsoil,LeafModel = 'PROSPECT-PRO')
        rdot<-data.foursail_pro[[1]]
        rsot<-data.foursail_pro[[2]]
        rfl.prosail<-ToolsRTM::Compute_BRF(rdot=rdot,rsot=rsot,tts=LUT[i,'tts'],data.light=ToolsRTM::dataSpec_PDB)
        sim.rfl[[i]]<-rfl.prosail
    
      } ##end paralle
    parallel::stopCluster(cl)
    sim.canopy<-do.call(rbind,sims)
    # remove box label
    sim.canopy <- sim.canopy[,-1] 
    if (is.null(method) | method == 'ggplot' | method != 'classical' ){
      
      sim.canopy.ggplot <-data.frame(Trait= LUT[,trait],ID=1:nrow(sim.canopy), sim.canopy)
      colnames(sim.canopy.ggplot)<-c('Trait','ID',paste('R',c(400:2499),sep='.'))
      sim.canopy.ggplot <-cbind(sim.canopy.ggplot[1], stack(sim.canopy.ggplot[3:dim(sim.canopy.ggplot)[2]]))
      sim.canopy.ggplot$ind =as.numeric(sub("R.", "", sim.canopy.ggplot$ind, fixed = TRUE))
      
      spectral_plot<-ggplot(sim.canopy.ggplot, aes(x=ind, y=values, group=Trait)) +theme_bw() +
        geom_line(aes(color=as.factor(Trait)),linetype = "dashed", size=0.8) +
        labs(color = trait, x='Wavelength',y='Reflectance')+  ggtitle("spectral signatures")+
        geom_point(aes(color=as.factor(Trait)), size=0.4) +
        theme(plot.title = element_text(hjust = 0.5, size=18))
      print(spectral_plot)
      
      to_export = list('LUT'=LUT,'Plot'=spectral_plot)
      
    } else if (method == 'classical') {
      
      # add the 'wavelength' and rotate the df
      # (i didn't find the actual wavelength values, but hey).
      sim.canopy <- cbind(1:ncol(sim.canopy), t(sim.canopy)) 
      plotspectra(sim.canopy, y=2:ncol(sim.canopy), cols = c(rainbow(ncol(sim.canopy)-2),'black'),
                  type='l', main=trait,ylab="reflectance", xlab = 'Wavelength')
      
      legend("topright", legend = seq(nmin, nmax, by=Interval),
             fill=c(rainbow(ncol(sim.canopy)-1),'black'),cex=0.8)
      
      to_export = LUT
    }
    
  } else if (model == 'PROSPECT') {
    
    message('PROPECT-PRO is processing ...')
    
    LUT<- data.frame(N,Cab,Car,Anth,Cbrown,EWT,LMA,alpha,Prot,CBC)
    nLUT = dim(LUT)[1]
    
    no_cores <- parallel::detectCores() - 2 
    cl <- parallel::makeCluster(no_cores)
    doParallel::registerDoParallel(cl)
    
    start_time <- Sys.time()
    sim.rfl<-list()
      ### simulations
      sims<-foreach(i=1:nLUT) %dopar% {
        rfl.prospect<-ToolsRTM::prospect_PRO(N=LUT[i,'N'],Cab=LUT[i,'Cab'],Car=LUT[i,'Car'],
                                           Anth=LUT[i,'Anth'],Cbrown=LUT[i,'Cbrown'],
                                           EWT= LUT[i,'EWT'],LUT[i,'LMA'],
                                           alpha=LUT[i,'alpha'],
                                           Prot= LUT[i,'Prot'],CBC=LUT[i,'CBC'])
        rho	 <- 	rfl.prospect[[2]]
      
        sim.rfl[[i]]<-rho
      
      } ##end paralle
    parallel::stopCluster(cl)
    

    sim.leaf<-do.call(rbind,sims)
    # remove box label
    sim.leaf <- sim.leaf[,-1] 
    
    if (is.null(method) | method == 'ggplot' | method != 'classical' ){
      
      sim.leaf.ggplot <-data.frame(Trait= LUT[,trait],ID=1:nrow(sim.leaf), sim.leaf)
      colnames(sim.leaf.ggplot)<-c('Trait','ID',paste('R',c(400:2499),sep='.'))
      sim.leaf.ggplot <-cbind(sim.leaf.ggplot[1], stack(sim.leaf.ggplot[3:dim(sim.leaf.ggplot)[2]]))
      sim.leaf.ggplot$ind =as.numeric(sub("R.", "", sim.leaf.ggplot$ind, fixed = TRUE))
    
      spectral_plot<-ggplot(sim.leaf.ggplot, aes(x=ind, y=values, group=Trait)) +theme_bw() +
          geom_line(aes(color=as.factor(Trait)),linetype = "dashed", size=0.8) +
          labs(color = trait, x='Wavelength',y='Reflectance')+  ggtitle("spectral signatures")+
          geom_point(aes(color=as.factor(Trait)), size=0.4) +
          theme(plot.title = element_text(hjust = 0.5, size=18))
      print(spectral_plot)
      
      to_export = list('LUT'=LUT,'Plot'=spectral_plot)
    } else if (method == 'classical') {
    
      # add the 'wavelength' and rotate the df
      # (i didn't find the actual wavelength values, but hey).
      sim.leaf <- cbind(1:ncol(sim.leaf), t(sim.leaf)) 
      plotspectra(sim.leaf, y=2:ncol(sim.leaf), cols = c(rainbow(ncol(sim.leaf)-2),'black'),
                  type='l', main=trait,ylab="reflectance", xlab = 'Wavelength')
      legend("topright", legend = seq(nmin, nmax, by=Interval),
             fill=c(rainbow(ncol(sim.leaf)-1),'black'),cex=0.8)
      
      to_export = LUT
    }
       
  } else if ( model == 'INFORM'){
      message('INFORM model is processing ...')
      LUT<- data.frame(N,Cab,Car,Anth,Cbrown,EWT,LMA,alpha,Prot,CBC,
                       LIDFa, LIDFb=0, TypeLidf=2, LAI=lai,
                       hspot=hotspot, tts, tto, psi,LAIu,sd,cd,h,skyl)
      nLUT = dim(LUT)[1]
      
      ## choose number of processors/cores
      no_cores <- parallel::detectCores() - 2 
      cl <- parallel::makeCluster(no_cores)
      doParallel::registerDoParallel(cl)
      
      start_time <- Sys.time()
      sim.rfl<-list()
      ### simulations
        sims<-foreach(i=1:nLUT) %dopar% {
          rfl.inform<-ToolsRTM::inform(inputLUT = LUT[i,],rsoil=rsoil,LeafModel = 'PROSPECT-PRO')
          sim.rfl[[i]]<-rfl.inform
          
        } ##end paralle
      
      parallel::stopCluster(cl)
      
      sim.canopy<-do.call(rbind,sims)
      # remove box label
      sim.canopy <- sim.canopy[,-1] 
      
      if (is.null(method) | method == 'ggplot' | method != 'classical' ){
        
        sim.canopy.ggplot <-data.frame(Trait= LUT[,trait],ID=1:nrow(sim.canopy), sim.canopy)
        colnames(sim.canopy.ggplot)<-c('Trait','ID',paste('R',c(400:2499),sep='.'))
        sim.canopy.ggplot <-cbind(sim.canopy.ggplot[1], stack(sim.canopy.ggplot[3:dim(sim.canopy.ggplot)[2]]))
        sim.canopy.ggplot$ind =as.numeric(sub("R.", "", sim.canopy.ggplot$ind, fixed = TRUE))
        
        spectral_plot<-ggplot(sim.canopy.ggplot, aes(x=ind, y=values, group=Trait)) +theme_bw() +
          geom_line(aes(color=as.factor(Trait)),linetype = "dashed", size=0.8) +
          labs(color = trait, x='Wavelength',y='Reflectance')+  ggtitle("spectral signatures")+
          geom_point(aes(color=as.factor(Trait)), size=0.4) +
          theme(plot.title = element_text(hjust = 0.5, size=18))
        print(spectral_plot)
        
        to_export = list('LUT'=LUT,'Plot'=spectral_plot)
        
      } else if (method == 'classical') {
        
        # add the 'wavelength' and rotate the df
        # (i didn't find the actual wavelength values, but hey).
        sim.canopy <- cbind(1:ncol(sim.canopy), t(sim.canopy)) 
        plotspectra(sim.canopy, y=2:ncol(sim.canopy), cols = c(rainbow(ncol(sim.canopy)-2),'black'),
                    type='l', main=trait,ylab="reflectance", xlab = 'Wavelength')
        
        legend("topright", legend = seq(nmin, nmax, by=Interval),
               fill=c(rainbow(ncol(sim.canopy)-1),'black'),cex=0.8)
        to_export = LUT
      }
      
      
     }
       
return(to_export)
  }
  

  
#' Plot one or more spectra from a data frame
#'
#' @param df data.frame. Table containing the spectra to plot, with wavelengths in one column and reflectance/radiance values in one or more other columns.
#' @param x integer or character. Column index or name in \code{df} to use as the x-axis (typically wavelength). Default 1.
#' @param y integer or character vector. Column index(es) or name(s) in \code{df} to plot as spectra on the y-axis. Default 2.
#' @param cols vector. Line colors, one per element of \code{y}. Defaults to \code{y} itself.
#' @param xlim numeric vector of length 2. X-axis plot limits. Defaults to the data range of column \code{x}.
#' @param ylim numeric vector of length 2. Y-axis plot limits. Defaults to the data range of columns \code{y}.
#' @param main character. Plot title.
#' @param xlab character. X-axis label.
#' @param ylab character. Y-axis label.
#' @param ... additional graphical parameters passed to the underlying \code{plot()} call.
#'
#' @return plot
#' @export
#'
plotspectra <- function(df, x=1, y=2, cols=y,
                          xlim=range(df[,x], na.rm = T),
                          ylim=range(df[,y], na.rm = T),
                          main="", xlab="", ylab="", ...){
    # setup plot frame
    par(mfrow=c(1,1),  mar = c(5,5,1.1,1),bg= "white", 
        font.main=1.5, cex.main=1.2,font.axis=2, cex.axis=1.2, las=1,
        font.lab=2, cex.lab=1.0)
    plot(NULL, 
         xlim=xlim, 
         ylim=ylim,
         main=main, xlab=xlab, ylab=ylab)
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
    grid(lwd = 2) 
    # plot all your y's against your x
    pb <- sapply(seq_along(y), function(i){
      points(df[,c(x, y[i])], col=cols[i], ...)
    })
  }
  
#' Gradient scale
#' Map a numeric vector to a color gradient
#'
#' @param x numeric vector. Values to map onto the color scale.
#' @param colors character vector. Colors defining the gradient, interpolated in order. Default \code{c("red","yellow","green")}.
#' @param colsteps integer. Number of discrete color steps in the gradient. Default 100.
#'
#' @return a scale
#' @export
#'
color.gradient <- function(x, colors=c("red","yellow","green"), colsteps=100) {
  
  return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
}
  
  
  

