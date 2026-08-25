
#' Simulate spectra with RT models
#'
#' @param inputLUT is a LUT with all the inputs needed for simulating a leaf or canopy model
#' @param psoil  soil factor, by default is set to 0.5
#' @param rtm.model model for generating the simulations ('PROSPECT-PRO','PROSAIL' or 'INFORM')
#'
#' @return a LUT with all the simulations and the plot
#' @export
#'
#' 
get_simulations<-function(inputLUT = NULL, psoil=0.5, rtm.model = 'PROSAIL'){
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
  

  if (is.null(LUT)){
    message('LUT is empty, please add a LUT table with the input for simulation')
    stop()
  }  
  
  if (is.null(rtm.model)){
    message('PROSPECT-PRO + fourSAIL will be used as  model')
    
  } else if (is.null(rtm.model) | rtm.model == 'PROSAIL'){
      message('coupled PROPECT + fourSAIL model is processing ...')
      
      if (dim(inputLUT)[2] < 18 ){
        message('please, check the LUT table, PROSAIL needs 18 inputs')
        stop()
      }
      #get leaf parameters from LUT
      
      N=inputLUT[,'N']; Cab=inputLUT[,'Cab']; Car=inputLUT[,'Car']; Anth=inputLUT[,'Anth']; Cbrown=inputLUT[,'Cbrown']
      EWT=inputLUT[,'EWT']; LMA=inputLUT[,'LMA']; alpha=inputLUT[,'alpha']
      Prot=inputLUT[,'Prot'];CBC=inputLUT[,'CBC']
      
      ## get canopy parameters from LUT
      LIDFa=inputLUT[,'LIDFa']; LIDFb=inputLUT[,'LIDFb']; TypeLidf=inputLUT[,'TypeLidf']; lai=inputLUT[,'LAI']
      hotspot=inputLUT[,'hspot']; tts=inputLUT[,'tts']; tto=inputLUT[,'tto']; psi=inputLUT[,'psi']
    
      nLUT = dim(inputLUT)[1]
      ## choose number of processors/cores
      no_cores <- max(1, parallel::detectCores() - 2)
      if (nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_"))) no_cores <- min(no_cores, 2L)
      cl <- parallel::makeCluster(no_cores)
      doParallel::registerDoParallel(cl)
      
      start_time <- Sys.time()
      sim.rfl<-list()
      ### simulations
        sims<-foreach::foreach(i=1:nLUT) %dopar% {
          data.foursail_pro<-ToolsRTM::foursail(inputLUT=inputLUT[i,],rsoil=rsoil,PROSPECTversion = 'PRO')
          rdot<-data.foursail_pro[[1]]
          rsot<-data.foursail_pro[[2]]
          rfl.prosail<-ToolsRTM::Compute_BRF(rdot=rdot,rsot=rsot,tts=inputLUT[i,'tts'],data.light=ToolsRTM::dataSpec_PDB)
          sim.rfl[[i]]<-rfl.prosail
          } ##end paralle
      parallel::stopCluster(cl)
      sim.canopy<-do.call(rbind,sims)
      # remove box label
      sim.canopy <- sim.canopy[,-1] 
      sim.canopy.to_export <-data.frame(inputLUT,ID=1:nrow(sim.canopy), sim.canopy)
      colnames(sim.canopy.to_export)<-c(names(inputLUT),'ID',paste('R',c(400:2499),sep='.'))
      
      names_withrfl <- names(sim.canopy.to_export[ , grepl( "R." , names(sim.canopy.to_export) ) ])

      sim.canopy.ggplot <-cbind(ID=sim.canopy.to_export[['ID']],group='sims', stack(sim.canopy.to_export[,names_withrfl]))
      sim.canopy.ggplot$ind =as.numeric(sub("R.", "", sim.canopy.ggplot$ind, fixed = TRUE))
      
      # create a new dataframe crop_means_Se
      group_means <- sim.canopy.ggplot |>  
        dplyr::group_by(ind, group) |> 
        dplyr::summarize(mean=mean(values), 
                  sd=sd(values), 
                  N=n(), 
                  quantile.25 = quantile(values, probs = 0.25),
                  quantile.75 = quantile(values, probs = 0.975),
                  se= sd / sqrt(N), 
                  upper_limit=mean + se, 
                  lower_limit=mean - se) 
      
      spectral_plot<-ggplot(group_means, aes(x=ind, y=mean,color=as.factor(group))) +theme_bw() +
        geom_line(aes(),linetype = "dashed", size=0.8) +
        geom_line(aes(y = quantile.25),linetype = "dashed") +
        geom_line(aes(y = quantile.75),linetype = "dashed") +
        labs(color=rtm.model, x='Wavelength',y='Reflectance')+  ggtitle("spectral simulations")+
        geom_point(aes(), size=0.4) +
        theme(plot.title = element_text(hjust = 0.5, size=18))
      print(spectral_plot)
      
      to_export = list('LUT'=sim.canopy.to_export,'Plot'=spectral_plot)
      
  } else if (rtm.model =='INFORM'){
    message('INFORM model is processing ...')
    
    if (dim(inputLUT)[2] < 24 ){
      message('please, check the LUT table, INFORM needs 24 inputs')
      stop()
    }
    #get leaf parameters from LUT
    
    N=inputLUT[,'N']; Cab=inputLUT[,'Cab']; Car=inputLUT[,'Car']; Anth=inputLUT[,'Anth']; Cbrown=inputLUT[,'Cbrown']
    EWT=inputLUT[,'EWT']; LMA=inputLUT[,'LMA']; alpha=inputLUT[,'alpha']
    Prot=inputLUT[,'Prot'];CBC=inputLUT[,'CBC']
    
    ## get canopy parameters from LUT
    LIDFa=inputLUT[,'LIDFa']; LIDFb=inputLUT[,'LIDFb']; TypeLidf=inputLUT[,'TypeLidf']; lai=inputLUT[,'LAI']
    hotspot=inputLUT[,'hspot']; tts=inputLUT[,'tts']; tto=inputLUT[,'tto']; psi=inputLUT[,'psi']
    
    ##  get canopy parameters from LUT needed for INFORM model
    
    laiu=inputLUT[,'LAIu']
    sd=inputLUT[,'sd']; cd=inputLUT[,'cd']; h=inputLUT[,'h']; psi=inputLUT[,'psi']
    skyl=inputLUT[,'skyl']
    
    nLUT = dim(inputLUT)[1]
    ## choose number of processors/cores
    no_cores <- parallel::detectCores() - 2 
    cl <- parallel::makeCluster(no_cores)
    doParallel::registerDoParallel(cl)
    
    start_time <- Sys.time()
    sim.rfl<-list()
    ### simulations
    sims<-foreach::foreach(i=1:nLUT) %dopar% {
      rfl.inform<-ToolsRTM::inform(inputLUT = inputLUT[i,],rsoil=rsoil,PROSPECTversion = 'PRO')
      sim.rfl[[i]]<-rfl.inform
    } ##end paralle
    parallel::stopCluster(cl)
    sim.canopy<-do.call(rbind,sims)
    # remove box label
    sim.canopy <- sim.canopy[,-1] 
    sim.canopy.to_export <-data.frame(inputLUT,ID=1:nrow(sim.canopy), sim.canopy)
    colnames(sim.canopy.to_export)<-c(names(inputLUT),'ID',paste('R',c(400:2499),sep='.'))
    
    names_withrfl <- names(sim.canopy.to_export[ , grepl( "R." , names(sim.canopy.to_export) ) ])
    sim.canopy.ggplot <-cbind(ID=sim.canopy.to_export[['ID']],group='sims', stack(sim.canopy.to_export[,names_withrfl]))
    sim.canopy.ggplot$ind =as.numeric(sub("R.", "", sim.canopy.ggplot$ind, fixed = TRUE))
    
    # create a new dataframe crop_means_Se
    group_means <- sim.canopy.ggplot |>  
      dplyr::group_by(ind,group) |> 
      dplyr::summarize(mean=mean(values), 
                sd=sd(values), 
                N=n(), 
                quantile.25 = quantile(values, probs = 0.25),
                quantile.75 = quantile(values, probs = 0.975),
    se= sd / sqrt(N), 
    upper_limit=mean + se, 
    lower_limit=mean - se ) 

    spectral_plot<-ggplot(group_means, aes(x=ind, y=mean,color=as.factor(group))) +theme_bw() +
      geom_line(aes(),linetype = "dashed", size=0.8) +
      geom_line(aes(y = quantile.25),linetype = "dashed") +
      geom_line(aes(y = quantile.75),linetype = "dashed") +
      labs(color=rtm.model, x='Wavelength',y='Reflectance')+  ggtitle("spectral simulations")+
      geom_point(aes(), size=0.4) +
      theme(plot.title = element_text(hjust = 0.5, size=18))
    print(spectral_plot)
    
    to_export = list('LUT'=sim.canopy.to_export,'Plot'=spectral_plot)
    
    } else if(rtm.model =='PROSPECT-PRO'){
      
      message(' PROSPECT-PRO  is processing ...')
      
      if (dim(inputLUT)[2] < 10 ){
        message('please, check the LUT table, PROSPECT-PRO needs 10 inputs')
        stop()
      }
      #get leaf parameters from LUT
    
      nLUT = dim(inputLUT)[1]
      ## choose number of processors/cores
      no_cores <- max(1, parallel::detectCores() - 2)
      if (nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_"))) no_cores <- min(no_cores, 2L)
      cl <- parallel::makeCluster(no_cores)
      doParallel::registerDoParallel(cl)
      
      start_time <- Sys.time()
      sim.rfl<-list()
      ### simulations
      sims<-foreach::foreach(i=1:nLUT) %dopar% {
        rfl.prospect<-ToolsRTM::prospect_PRO(N=inputLUT[i,'N'],Cab=inputLUT[i,'Cab'],Car=inputLUT[i,'Car'],
                                             Anth=inputLUT[i,'Anth'],Cbrown=inputLUT[i,'Cbrown'],
                                             EWT= inputLUT[i,'EWT'],inputLUT[i,'LMA'],
                                             alpha=inputLUT[i,'alpha'],
                                             Prot= inputLUT[i,'Prot'],CBC=inputLUT[i,'CBC'])
        rho	 <- 	rfl.prospect[[2]]
        
        sim.rfl[[i]]<-rho
      } ##end paralle
      parallel::stopCluster(cl)
      sim.leaf<-do.call(rbind,sims)
      # remove box label
      sim.leaf <- sim.leaf[,-1] 
      sim.leaf.to_export <-data.frame(inputLUT,ID=1:nrow(sim.leaf), sim.leaf)
      colnames(sim.leaf.to_export)<-c(names(inputLUT),'ID',paste('R',c(400:2499),sep='.'))
      
      names_withrfl <- names(sim.leaf.to_export[ , grepl( "R." , names(sim.leaf.to_export) ) ])
      
      sim.leaf.ggplot <-cbind(ID=sim.leaf.to_export[['ID']],group='sims', stack(sim.leaf.to_export[,names_withrfl]))
      sim.leaf.ggplot$ind =as.numeric(sub("R.", "", sim.leaf.ggplot$ind, fixed = TRUE))
      
      # create a new dataframe crop_means_Se
      group_means <- sim.leaf.ggplot |>  
        dplyr::group_by(ind,group) |> 
        dplyr::summarize(mean=mean(values), 
                  sd=sd(values), 
                  N=n(), 
                  quantile.25 = quantile(values, probs = 0.25),
                  quantile.75 = quantile(values, probs = 0.975),
      se= sd / sqrt(N), 
      upper_limit=mean + se, 
      lower_limit=mean - se ) 

  spectral_plot<-ggplot(group_means, aes(x=ind, y=mean,color=as.factor(group))) +theme_bw() +
    geom_line(aes(),linetype = "dashed", size=0.8) +
    geom_line(aes(y = quantile.25),linetype = "dashed") +
    geom_line(aes(y = quantile.75),linetype = "dashed") +
    labs(color=rtm.model, x='Wavelength',y='Reflectance')+  ggtitle("spectral simulations")+
    geom_point(aes(), size=0.4) +
    theme(plot.title = element_text(hjust = 0.5, size=18))
  print(spectral_plot)
  
  to_export = list('LUT'=sim.leaf.to_export,'Plot'=spectral_plot)
}
  
  return(to_export)
}
  
