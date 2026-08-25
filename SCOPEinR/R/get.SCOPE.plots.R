


#' \code{get.SCOPE.plots} get simulations based on SCOPE model


#' @param path.files  folder with main tables from  the SCOPE simulation.
#' @param plant.trait  plant trait
#' @param get.plots  type plot for comparing, options are: 'fluorescence', 'reflectance', 'radiance'
#' @author Carlos Camino (Ported version into R)
#'
#' @examples
#' \dontrun{
#' get.SCOPE.plots(path.files = "outs", plant.trait = "Cab", get.plots = "reflectance")
#' }
#' @export


get.SCOPE.plots <-function(path.files, plant.trait, get.plots)  {

  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ### 0.1 Check the get.plots ---
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  if (missing(get.plots)){

    type.plots = 'reflectance'

  } else{

    type.plots = get.plots

  }


  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  #: 1. Load main variables
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  if (missing(path.files)) {
    stop('Please specify the path with the output files')
  }

  path_to.plot <- paste(path.files,'/Plots.checks/',sep='')

  if (!dir.exists(path_to.plot)) {

      dir.create(path_to.plot)
  }

  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  #:: find files ----
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  ## LUT table
  files.lut <- list.files(paste(path.files,'/Parameters',sep = ''), pattern = "LUT", full.names = TRUE)

  if (type.plots == 'reflectance') {

    # Get the list of files that match the pattern
    files.rdd <- list.files(path.files, pattern = "rdd", full.names = TRUE)
    files.rdo <- list.files(path.files, pattern = "rdo", full.names = TRUE)
    files.rsd <- list.files(path.files, pattern = "rsd", full.names = TRUE)
    files.rso <- list.files(path.files, pattern = "rso", full.names = TRUE)
    files.rfl <- list.files(path.files, pattern = "refl", full.names = TRUE)
    files.all <- c(files.rdd,files.rdo,files.rsd,files.rso,files.rfl)

    name.plot <- '1-reflectance_'

  } else if (type.plots == 'fluorescence') {

    files.all <- list.files(path.files, pattern = "fluorescence", full.names = F)

    files.all <- setdiff(files.all,c("fluorescence_scalar_units.csv", "fluorescence_scalar.csv"))
    files.all <-paste(path.files,'/',files.all,sep='')

    name.plot <- '6-fluorescence_'

  } else if ( (type.plots == 'radiance') ){

    files.all <- list.files(path.files, pattern = "Lo_", full.names = T)

    name.plot <- '3-Radiance_'
  }

  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # Read each file as a CSV
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  data <- lapply(files.all, read.csv)
  lut <- lapply(files.lut, read.csv)

  if (missing(plant.trait)){

    message('Please select an input form LUT: a vector with more inputs is available.')
    stop()

  }

  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  #:: make plots ----
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  for (j in c(1:length(plant.trait))) {

    i.trait = plant.trait[j]
    print(i.trait)


    min.plant.trait <- min(lut[[1]][,i.trait], na.rm=T)
    max.plant.trait <- max(lut[[1]][,i.trait], na.rm=T)

    if (min.plant.trait == max.plant.trait){
      stop(paste('There is not variation in range for ',i.trait,sep=''))
    }



    # Access the data from each file

    for (i in seq_along(data)) {
      filenames <- basename(files.all[i])
      # Remove the extension ".csv" from each file name
      filename_ <- gsub("\\.csv$", "", filenames)
      df <- data[[i]]
      ndim_ <-dim(df)[1] - 0

      if (ndim_ < 10){
        stop(' A LUT with more than 10 variations in the inputs is needed')
      }

      df <- data.frame(trait = lut[[1]][,i.trait],
                       nsim = c(1:dim(df)[1]),
                       df)


      data.spectra<-SCOPEinR::define.bands()

      df.plot <- df %>%
        pivot_longer(cols = starts_with("X"),
                     names_to = "Wavelength",
                     values_to = "Spectrum")

      if (type.plots == 'fluorescence') {

        df.plot$Wavelength <- rep(data.spectra$wlF,ndim_)
        XLimit <- 640
        Ylimit <- 850
        Ylabel <- 'Fluorescence'

      } else if (type.plots == 'reflectance')  {
        ### 400- thermal range (2162 waves)
        df.plot$Wavelength <- rep(data.spectra$wlS[1:2001],ndim_)
        XLimit <- 400
        Ylimit <- 2499
        Ylabel <- 'Reflectance'

      } else {

        ### 400- thermal range (2162 waves)
        df.plot$Wavelength <- rep(data.spectra$wlS[1:2001],ndim_)
        XLimit <- 400
        Ylimit <- 2499
        Ylabel <- 'Radiance'

      }


      # Define the minimum and maximum values of trait for grouping
      trait.min <- round(min(df.plot[['trait']],na.rm=T),2)
      trait.median <- round(median(df.plot[['trait']],na.rm=T),2)
      trait.max <- round(max(df.plot[['trait']],na.rm=T),2)

      # Create breaks based on the min and max trait values
      num_groups <- 8
      breaks <- round(seq(trait.min, trait.max, length.out = num_groups+1),2)
      # Assign trait groups based on breaks
      df.plot[['group.trait']] <- cut(df.plot[['trait']], breaks = breaks,labels = F)#, paste("trait.", c(1:num_groups),sep=''))


      # Calculate average spectra for each group
      average.spectra <- df.plot %>%
        group_by(group.trait, Wavelength) %>%
        summarize(Average_Spectrum = mean(Spectrum, na.rm = TRUE)) %>%
        ungroup()

      # Create a color palette with custom colors for the groups
      custom_colors <- c("#FF0000", "#00FF00", "#0000FF", "#FF00FF", "#FFFF00")
      #custom_colors <- c("#0000FF","#FFFFFF", "#FF0000")

      # Plot the average spectra
      plot.trait <-ggplot(average.spectra, aes(x = Wavelength, y = Average_Spectrum, color =group.trait)) +
        geom_line(linewidth = 0.5) +

        scale_color_gradientn(name =i.trait,
                              colors = terrain.colors(10), #terrain.colors rainbow
                              guide = "colourbar",
                              space = "Lab",
                              breaks = breaks[-length(breaks)],
                              labels = NULL) +
        theme_bw() +   xlim(XLimit, Ylimit) +
        theme(legend.position = "none") + ylab(Ylabel ) + xlab('Wavelength') +
        guides(color = guide_legend(title =i.trait ))

      # Display the plot
      ggsave(paste(path_to.plot,name.plot,filename_,'_by_',i.trait,'.png',sep=''),
             plot=plot.trait,
             width = 10, height = 10,  dpi = 300,units = "cm")

    } #end for plot
    ###############################################################
    ###############################################################
  } ## end for i.trait
  ###############################################################



}












