
rm(list= ls())

# 1. load the main libraries  -----

library(ToolsRTM)
required_packages <- c("shiny", "shinythemes", "ggplot2", "dplyr", 
                       "doParallel",'foreach','DT','parallel')

# Check for missing packages and install them if necessary
missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}
# Load the libraries
lapply(required_packages, library, character.only = TRUE)

# Define UI

ui <- navbarPage("Online reflectance simulator",theme = shinytheme("flatly"),

                 tabPanel(title = "BRDF",
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    # Sidebar panel for inputs ----

    # Sidebar panel for inputs ----
    sidebarPanel(width = 6,

                 fluidRow(h4("Leaf parameters"),

                          column(width = 4,

                                 sliderInput(inputId = "Cab",
                                             label = "Chlorophyll content",
                                             min = 0,
                                             max = 80,
                                             value = 40),
                                 sliderInput(inputId = "Cbrown",
                                             label = "Brown pigments",
                                             min = 0,
                                             max = 1,
                                             step = 0.05,
                                             value = 0.0),


                          ),

                          column(width = 4,
                                 sliderInput(inputId = "Car",
                                             label = "Carotenoid content",
                                             min = 0,
                                             max = 20,
                                             value = 8),
                                 sliderInput(inputId = "EWT",
                                             label = "Equivalent water thickness",
                                             min = 0.0001,
                                             max = 0.05,
                                             value = 0.01),


                          ),

                          column(width = 4,
                                 sliderInput(inputId = "Anth",
                                             label = "Anthocyanin Content",
                                             min = 0,
                                             max = 7,
                                             step = 0.2,
                                             value = 2),
                                 sliderInput(inputId = "CBC",
                                             label = "Carbon-based constituent",
                                             min =  0.001,
                                             max = 0.04,
                                             value = 0.03),

                          ),
                          column(width = 4,
                                 sliderInput(inputId = "N",
                                             label = "Structure parameter N",
                                             min = 1,
                                             max = 3,
                                             step = 0.2,
                                             value = 1.5),


                          ),


                          column(width = 4,

                                 sliderInput(inputId = "Prot",
                                             label = "Leaf protein content",
                                             min = 0.0001,
                                             max = 0.005,
                                             value = 0.0001)
                          )),

                 fluidRow(h4("Canopy parameters"),

                          column(width = 4,

                                 sliderInput(inputId = "LAI",
                                             label = "Leaf Area Index",
                                             min = 0.001,
                                             max = 10,
                                             step = 0.1,
                                             value = 4)),

                          column(width = 4,

                                 sliderInput(inputId = "hspot",
                                             label = "Hotspot parameter",
                                             min = 0,
                                             max = 1,
                                             value = 0.01)),

                          column(width = 4,

                                 selectInput(inputId = "TypeLIDF",
                                             label = "Type LIDF",
                                             choices = list(Planophile = "plano",
                                                            Erectophile = "erecto",
                                                            Plagiophile = "plagio",
                                                            Extremophile = "extremo",
                                                            Uniform = "uniform",
                                                            Spherical = "sph"),
                                             selected = "sph"))
                 ),

                 fluidRow(h4("Soil and angles parameters"),

                          column(width = 4,

                                 sliderInput(inputId = "psoil",
                                             label = "Soil brightness",
                                             min = 0,
                                             max = 1,
                                             value = 0.5)),
                          column(width = 4,

                                 sliderInput(inputId = "wl",
                                             label = "Wavelength [nm]",
                                             min = 400,
                                             max = 2500,
                                             value = 550)),

                          column(width = 4,

                                 radioButtons(inputId = "angle_step",
                                              label = "Angular resolution",
                                              choices = list("2.5" = 2.5,
                                                             "5" = 5,
                                                             "10" = 10,
                                                             "20" = 20),
                                              selected = 20)),
                 ),

                 fluidRow(h4("Viewing angle parameters"),

                          column(width = 4,

                                 sliderInput(inputId = "tts",
                                             label = "Solar zenith angle",
                                             min = 0,
                                             max = 90,
                                             value = 30)),
                          column(width = 4,

                                 sliderInput(inputId = "tto",
                                             label = "Observer zenith angle",
                                             min = 0,
                                             max = 90,
                                             value = 10)),

                          column(width = 4,

                                 sliderInput(inputId = "psi",
                                             label = "Relative azimuth angle",
                                             min = 0,
                                             max = 180,
                                             value = 0)),

                          ),

                 p('Note: Adjustments to the observer zenith angle (tto) and PSI values do not impact on the simulations.
                   All simulations cover angles from 0 degrees (nadir) to 90 degrees in both forward (positive)
                   and backward (negative) directions, based on the provided angular resolution.'),


                 # Add options to visualize reflectance and transmittance at leaf level
                 fluidRow(h4("Leaf Optical Properties"),

                          column(width = 4,
                                 checkboxInput(inputId = "show_reflectance",
                                               label = "Show Reflectance at Leaf Level",
                                               value = FALSE)),


                          column(width = 4,
                                 checkboxInput(inputId = "show_transmittance",
                                               label = "Show Transmittance at Leaf Level",
                                               value = FALSE)),
                          column(width = 4,
                                 checkboxInput(inputId = "show_diff",
                                               label = "Show spectral differences",
                                               value = FALSE))
                 )
    ),

    # Main panel for displaying outputs ----
    mainPanel(width = 6,

              # Create a tabset panel
              tabsetPanel(id = "tabs",
                tabPanel("Simulations",
                         h4("Simulated Spectrum"),
                         div(style = "overflow-y: auto; height: 400px;",  # Set a fixed height and enable vertical scrolling
                             plotOutput(outputId = "simulations_withAngles")
                         ),
                         DT::dataTableOutput("dataset_LUT"),
                         div(style = "overflow-y: auto; height: 400px;",  # Combined plot (reflectance and transmittance)
                             plotOutput(outputId = "refl_trans_plot"))
                ),
                tabPanel("Spectral Differences",
                         div(style = "overflow-y: auto; height: 800px;",  # Set a fixed height for the diff plot
                             plotOutput(outputId = "diff_plot",height = "600px") )

                ),
                tabPanel("BRDF",
                         div(style = "overflow-y: auto; height: 1200px;",  # Set a fixed height for the diff plot
                             plotOutput(outputId = "simulations_brdf",height = "600px" ))
                ),
                tabPanel("about the model",

                         # Introduction
                         h4("Introduction"),
                         p(style = "text-align: justify;", "This module is designed to explore various leaf inclination distributions
                                  observed in crop and forest canopies, including erectophile, planophile, spherical,
                                  plagiophile, and extremophile distributions. By using the PROSPECT-PRO and fourSAILH models,
                                  we analyze how viewing angles (such as sun zenith angle (tts) and observer zenith angle (tto)) and relative azimuth angle (psi)
                                    affect reflectance across wavelengths (400-2500 nm)."),

                         p(style = "text-align: justify;", "The module allows for the examination of different angular resolutions (2.5°, 5°, and 10°),
                                    helping users understand the influence of these angles on the bidirectional reflectance distribution
                                    function (BRDF) properties. It also provides insights into the variation of reflectance when observed
                                    from different view zenith angles (VZA) along the principal plane (the plane containing the sun, the observer, and the target)
                                    and the cross-principal plane, demonstrating how reflectance changes with viewing geometry."),
                         p(style = "text-align: justify;", "the 'polar plot' (in BRDF tab), shows the view zenith angle (VZA) runs from the center (nadir, 0°) to the outer circle (90°),
                                    while the relative azimuth angle (psi) changes along the circumference from 0° to 360°, with the sun positioned
                                    at psi equal to 0°. This visualization provides an intuitive way to analyze the BRDF and its
                                    spherical properties, helping users see how reflectance varies with different view geometries."),

                         # PROSPECT-PRO Model Table
                         h4("PROSPECT-PRO inputs"),

                         tableOutput("prospect_pro_table_brf"),

                         # fourSAILH Model Table
                         h4("fourSAILH inputs"),

                         tableOutput("foursailh_table_brf"),

                         # Leaf Inclination Distribution Table
                         h4("Leaf Inclination Distribution Function (LIDF) Types"),

                         p(style = "text-align: justify;",
                           "The table below shows different leaf inclination distributions observed in
                                    plant canopies and their corresponding LIDFa and LIDFb values:"),

                         tableOutput("lidf_table")
                )# end tabPanel
              )

    )
  )
 )
)
# Define server logic required to draw a histogram ----
server <- function(input, output,session) {

  # Histogram of the Old Faithful Geyser Data ----
  # with requested number of bins
  # This expression that generates a histogram is wrapped in a call
  # to renderPlot to indicate that:
  #
  # 1. It is "reactive" and therefore should be automatically
  #    re-executed when inputs (input$bins) change
  # 2. Its output type is a plot


  # Define the reactive expression for the leaf parameters ---------------------------------------------------
  parameters_model <- reactive({
    params <- c(
      N = as.numeric(input$N),
      Cab = as.numeric(input$Cab),Car = as.numeric(input$Car),Anth = as.numeric(input$Anth),
      Cbrown = as.numeric(input$Cbrown),EWT = as.numeric(input$EWT),
      LMA = 0,                # Assuming LMA is constant, no need for conversion
      alpha = 40,             # Same with alpha
      Prot = as.numeric(input$Prot),CBC = as.numeric(input$CBC),
      LAI = as.numeric(input$LAI),
      TypeLidf = 1,           # Assuming 1 for considering LIDFb=0. and LIDFa in degrees.
      hspot = as.numeric(input$hspot),
      tts = as.numeric(input$tts),tto = as.numeric(input$tto),psi = as.numeric(input$psi),
      psoil = as.numeric(input$psoil))

    return(params)
  })

  # Reactive expression for the LUT data ---------------------------------------------------
  lut_data.sim <- reactive({
    # Get the selected leaf and canopy parameters from the reactive `params` expression
    lut.to_sim <- data.frame(t(as.data.frame(parameters_model())))
    row.names(lut.to_sim) <- NULL


    # Call the Simulations
    data <- ToolsRTM::dataSpec_PDB
    Rsoil.dry  <- data[,11]  # rsoil1 = dry soil
    Rsoil.wet <- data[,12]   # rsoil2 = wet soil
    print(lut.to_sim$psoil)
    psoil <- lut.to_sim[1,'psoil']

    # soil factor (psoil=0: wet soil / psoil=1: dry soil)
    rsoil <- psoil * Rsoil.dry + (1 - psoil) * Rsoil.wet




    return(list(lut.to_sim, rsoil))
  })


  # Observe the checkbox input to switch to the "Spectral Differences" tab
  observeEvent(input$show_diff, {
    if (input$show_diff) {
      updateTabsetPanel(session, "tabs", selected = "Spectral Differences")
    }
  })

  # Observe the checkbox input to switch to the "Spectral Differences" tab
  observeEvent(input$show_reflectance, {
    if (input$show_reflectance) {
      updateTabsetPanel(session, "tabs", selected = "Simulations")
    }
  })

  # Observe the checkbox input to switch to the "Spectral Differences" tab
  observeEvent(input$show_transmittance, {
    if (input$show_transmittance) {
      updateTabsetPanel(session, "tabs", selected = "Simulations")
    }
  })

  # Reactive expression for reflectance data ---------------------------------------------------
  reflectance_data_prosail <- reactive({
    # Show the modal window
    shinybusy::show_modal_spinner()


    df.LUT <- lut_data.sim()[[1]]
    row.names(df.LUT) <- NULL
    print(df.LUT)

    rsoil <- lut_data.sim()[[2]]

    # Update LIDF parameters based on input
    df.LUT$LIDFa <- switch(input$TypeLIDF,
                           plano = 1,
                           erecto = -1,
                           plagio = 0,
                           extremo = 0,
                           uniform = 0,
                           sph = -0.35)

    df.LUT$LIDFb <- switch(input$TypeLIDF,
                           plano = 0,
                           erecto = 0,
                           plagio = -1,
                           extremo = 1,
                           uniform = 0,
                           sph = -0.15)


    print(df.LUT)

    angle_step <- as.numeric(input$angle_step)

    angles <- expand.grid(tto = seq(0, 90 - angle_step, angle_step),
                          psi = seq(0, 180, angle_step))
    wl.selected = input$wl

    # Assuming LUT_ is your DataFrame
    df.LUT <- df.LUT[, !(names(df.LUT) %in% c("tto", "psi"))]
    df.LUT <-data.frame(df.LUT,angles)
    dim(df.LUT)
    wavelength <- seq(400,2500,1)

    # choose number of processors/cores
    no_cores <- parallel::detectCores() - 2
    cl <- parallel::makeCluster(no_cores)
    doParallel::registerDoParallel(cl)
    rho<-list()
    tau<-list()
    # Export df.LUT and other required objects to the workers
    parallel::clusterExport(cl, varlist = c("df.LUT"), envir = environment())

    sims<-foreach::foreach(i=1:dim(df.LUT)[1]) %dopar% {

      sim_leaf_values <- ToolsRTM::prospect_PRO(N=df.LUT[i,'N'],Cab=df.LUT[i,'Cab'], Car=df.LUT[i,'Car'], Anth=df.LUT[i,'Anth'],
                                                Cbrown=df.LUT[i,'Cbrown'],
                                                EWT=df.LUT[i,'EWT'],LMA=df.LUT[i,'LMA'],alpha=40,
                                                Prot=df.LUT[i,'Prot'],CBC=df.LUT[i,'CBC'])
      wavelength <- sim_leaf_values[[1]]
      rho	 <- 	sim_leaf_values[[2]]
      tau	 <- 	sim_leaf_values[[3]]

      data <- ToolsRTM::dataSpec_PDB
      Rsoil.dry  <- data[,11]  # rsoil1 = dry soil
      Rsoil.wet <- data[,12]   #
      psoil = df.LUT[i,'psoil']
      rsoil<- c(psoil*Rsoil.dry+(1-psoil)*Rsoil.wet)

      data.foursail<- ToolsRTM::foursail(inputLUT=df.LUT[i,],rsoil=rsoil,LeafModel = 'PROSPECT-PRO')
      rdot<-data.foursail[[1]]
      rsot<-data.foursail[[2]]
      rfl.brf<-ToolsRTM::Compute_BRF(rdot=rdot,rsot=rsot,tts=df.LUT[i,'tts'],data.light=ToolsRTM::dataSpec_PDB)
      list(wavelength = wavelength, rho = rho, tau = tau, rfl = rfl.brf)
    } ##end paralle

    # Close the cluster after parallel processing is complete
    parallel::stopCluster(cl)

    #sim.canopy<-do.call(rbind,sims)
    # Combine the results and calculate the averages
    rho_matrix <- do.call(cbind, lapply(sims, function(x) x$rho))
    tau_matrix <- do.call(cbind, lapply(sims, function(x) x$tau))
    sim.canopy <- do.call(cbind, lapply(sims, function(x) x$rfl))

    # Calculate the average across all simulations for each wavelength
    avg_rho <- rowMeans(rho_matrix)
    avg_tau <- rowMeans(tau_matrix)

    # Plot reflectance/ transmittance at leaf level  ---------------------------------------------------
    output$refl_trans_plot <- renderPlot({
      # Check if either reflectance or transmittance is selected
      if (input$show_reflectance || input$show_transmittance) {

        # Initialize an empty dataframe for combined plotting
        df_combined <- data.frame()

        # Add reflectance data if selected
        if (input$show_reflectance) {
          df_combined <- rbind(
            df_combined,
            data.frame(wavelength = wavelength, value = avg_rho, variable = "Reflectance")
          )
        }

        # Add transmittance data if selected
        if (input$show_transmittance) {
          df_combined <- rbind(
            df_combined,
            data.frame(wavelength = wavelength, value = avg_tau, variable = "Transmittance")
          )
        }
        # Show a notification when the simulation data is processed
        showNotification("Leaf-level simulation successfully completed.", type = "message")


        # Plot both reflectance and transmittance in a single plot
        ggplot(df_combined, aes(x = wavelength, y = value, color = variable)) +
          geom_line(linewidth = 1) +
          facet_grid(variable ~ ., scales = "free_y") +  # Create facets for reflectance and transmittance
          scale_color_manual(values = c("Reflectance" = "forestgreen", "Transmittance" = "navyblue")) +
          labs(x = "Wavelength (nm)", y = NULL) +  # No common y-axis label
          theme_bw() +
          theme(
            legend.position = 'none',
            strip.text = element_text(face = "bold", size = 14),  # Bold facet labels
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            axis.title = element_text(face = "bold", size = 14),
            axis.text.y = element_text(hjust = 0.5, size = 12, face = "bold"),
            axis.text.x = element_text(hjust = 0.5, size = 12, face = "bold"),
            panel.grid.major = element_blank(),  # Optional: Remove grid lines
            panel.grid.minor = element_blank())
      }
    })

    # Plot spectral differences between canopy and leaf  ---------------------------------------------------


    # Notify the user that the simulation may take time
    showNotification("Starting forward simulation, this may take time depending on the matrix generated by the angular resolution.", type = "warning")

    # choose number of processors/cores
    no_cores <- parallel::detectCores() - 2
    cl <- parallel::makeCluster(no_cores)
    doParallel::registerDoParallel(cl)
    sim.rfl<-list()

    sims<-foreach::foreach(i=1:dim(df.LUT)[1]) %dopar% {
      data <- ToolsRTM::dataSpec_PDB
      Rsoil.dry  <- data[,11]  # rsoil1 = dry soil
      Rsoil.wet <- data[,12]   #
      psoil = df.LUT[i,'psoil']
      rsoil<- c(psoil*Rsoil.dry+(1-psoil)*Rsoil.wet)

      data.foursail<- ToolsRTM::foursail(inputLUT=df.LUT[i,],rsoil=rsoil,LeafModel = 'PROSPECT-PRO')
      rdot<-data.foursail[[1]]
      rsot<-data.foursail[[2]]
      rfl.prosail<-ToolsRTM::Compute_BRF(rdot=rdot,rsot=rsot,tts=df.LUT[i,'tts'],data.light=ToolsRTM::dataSpec_PDB)
      sim.rfl[[i]]<-rfl.prosail

    } ##end paralle

    # Close the cluster after parallel processing is complete
    parallel::stopCluster(cl)

    sim.canopy<-do.call(rbind,sims)

    # Calculate the average across all simulations for each wavelength
    reflectance_values <- colMeans(sim.canopy)

    shinybusy::remove_modal_spinner()


    output$diff_plot <- renderPlot({
      if (input$show_diff) {  # Add checkbox input for showing differences

        # Create the dataframe for differences
        df.diff<- data.frame(
          wavelength = wavelength,
          refl.canopy = reflectance_values,
          refl.leaf = avg_rho,
          diff = reflectance_values - avg_rho,
          diff_sqrt = (reflectance_values - avg_rho)^2,
          percent_diff = (reflectance_values - avg_rho) / avg_rho * 100)

        # Reshape the data for plotting with facets
        df.melted <- reshape2::melt(df.diff, id.vars = "wavelength",
                                    measure.vars = c("diff", 'diff_sqrt',"percent_diff"),
                                    variable.name = "variable",
                                    value.name = "value")





        # Define custom labels for facets
        custom_labels <- c(
          "diff" = "Absolute Difference",
          'diff_sqrt' = 'Squared Difference',
          "percent_diff" = "Difference in %"
        )

        # Plot using facets
        plot.dif <- ggplot(df.melted, aes(x = wavelength, y = value, color = variable)) +
          geom_line(linewidth = 1) +
          facet_grid(variable ~ ., scales = "free_y", labeller = as_labeller(custom_labels)) +  # Create facets with custom labels
          scale_color_manual(values = c("diff" = "gray6", "percent_diff" = "dodgerblue4", 'diff_sqrt' = 'darkolivegreen4')) +
          labs(x = "Wavelength (nm)",
               y = NULL,
               title = "Difference Between Canopy and Leaf Reflectance") +
          theme_bw() +
          theme(
            legend.position = 'none',
            strip.text = element_text(face = "bold", size = 14),  # Bold facet labels
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            axis.title = element_text(face = "bold", size = 14),
            axis.text.y = element_text(hjust = 0.5, size = 12, face = "bold"),
            axis.text.x = element_text(hjust = 0.5, size = 12, face = "bold"),
            panel.grid.major = element_blank(),  # Optional: Remove grid lines
            panel.grid.minor = element_blank()
          )


        # Show notification when the plot is generated
        showNotification("Difference between canopy and leaf reflectance successfully generated.", type = "message")
        #print(plot.dif)
        print(plot.dif)
      }
    })


    return(list(sim.canopy,wavelength,angles,df.LUT))



  })


  # Plot reflectance/ transmittance at leaf level  ---------------------------------------------------
  output$simulations_withAngles <- renderPlot({
    wavelength <-reflectance_data_prosail()[[2]]
    sim.canopy <- reflectance_data_prosail()[[1]]
    df <- data.frame(sim.canopy)

    # Reshape the data to long format
    df_long <- tidyr::gather(df, key = "band", value = "value")
    # Make 'band' an ordered factor with desired order
    df_long$band <- factor(df_long$band, levels = paste0("X", 1:ncol(df)))
    # Calculate average, 25th percentile, and 50th percentile for each band
    summary_stats <- df_long %>%
      group_by(band) %>%
      summarise(
        average = mean(value),
        median = median(value),
        percentile_25 = quantile(value, 0.25),
        percentile_50 = quantile(value, 0.50),
        percentile_75 = quantile(value, 0.75)
      )
    summary_stats$band <-wavelength
    # Include the number of rows in the notification message
    num_rows <- nrow(df)
    showNotification(paste("Forward simulation done successfully for", num_rows, "simulations"), type = "message")

    # Plot using ggplot2
    plot.sim <-ggplot(summary_stats, aes(x = band)) +
      geom_line(aes(y = average), color = "black", linewidth = 0.6) +
      geom_line(aes(y = median),  linetype = "dashed", color = "black", linewidth = 0.6) +
      geom_ribbon(aes(ymin = percentile_25, ymax = percentile_75), linetype = "dashed",fill = "black", alpha = 0.3) +
      labs(
        title = "",
        x = "wavelength (nm)",
        y = "Reflectance"
      ) +
      theme_bw()  +  theme(
        text = element_text(size = 14, face='bold'),  # Increase the text size
        axis.title = element_text(size = 16, face = "bold"),  # Make axis titles bold
        axis.text = element_text(face = "bold"),  # Make axis numbers bold
        plot.title = element_text(face = "bold"),  # Make plot title bold
        plot.subtitle = element_text(size = 14, face='bold')  # Adjust subtitle size
      )

    print(plot.sim)

    showNotification("LUT with simulations done successfully.", type = "message")
  })


  # Plot reflectance/ transmittance at leaf level  ---------------------------------------------------
  output$simulations_brdf <- renderPlot({
    sim.canopy <- reflectance_data_prosail()[[1]]
    df <- data.frame(sim.canopy)
    angles <- reflectance_data_prosail()[[3]]
    print(head(angles))
    angle_step <- as.numeric(input$angle_step)
    print(angle_step)
    wl.selected <- input$wl - 399
    wl.i <-paste('X',wl.selected,sep = '')
    print(wl.selected)
    # Ensure that the selected wavelength is within the range of your data (400 to 2500 nm)
    if (input$wl >= 400 && input$wl <= 2500) {
      # Find the corresponding column in the df (assuming column names are actual wavelengths)
      col_index <- which(names(df) == wl.i)
      print(col_index)
    }

    # Generate data for plotting
    print(head(df[1:3,1:10]))
    dat <- data.frame(BRDF = df[, col_index]) %>%
      cbind(angles) %>%
      mutate(psi_min = pmax(0, psi - angle_step / 2),
             psi_max = pmin(180, psi + angle_step / 2),
             tto_min = tto - angle_step / 2,
             tto_max = tto + angle_step / 2)
    print(head(dat))

    dat_fullpsi <- rbind(dat,
                         dat %>%
                           mutate(psi_min = 360 - psi_min,
                                  psi_max = 360 - psi_max))


    brdf <- ggplot(dat_fullpsi, aes(xmin = psi_min, xmax = psi_max, ymin = tto_min, ymax = tto_max, fill = BRDF)) +
      geom_rect() +
      annotate(geom = "segment", x = 0, y = 90, xend = 0, yend = 0) +
      annotate(geom = "segment", x = 180, y = 0, xend = 180, yend = 90) +
      annotate(geom = "segment", x = 90, y = 0, xend = 90, yend = 90, linetype = "dashed") +
      annotate(geom = "segment", x = 270, y = 0, xend = 270, yend = 90, linetype = "dashed") +
      coord_polar(theta = "x") +
      theme_bw() +
      xlab("PSI") +
      ylab("VZA") +
      scale_x_continuous(limits = c(0, 360), breaks = seq(0, 360, by = 45), expand = c(0, 0)) +
      scale_fill_viridis_c() +
       theme( text = element_text(size = 14, face='bold'),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 14, face='bold') )


    # principal plane plots
    pp_dat <- dat %>%
      dplyr::filter(psi == 0 | psi == 180) %>%
      mutate(tto = ifelse(psi > 0, tto, -tto))

    pp <- ggplot(pp_dat, aes(x = tto, y = BRDF)) +
      geom_line() +
      xlab("VZA") +
      ylab("BRDF") +
      theme_bw() +
      ggtitle("Principal plane") +
      scale_x_continuous(limits = c(-90, 90), breaks = seq(-90, 90, by = 45), expand = c(0, 0))+
      theme( text = element_text(size = 14, face='bold'),
             axis.title = element_text(size = 16, face = "bold"),
             axis.text = element_text(face = "bold"),
             plot.title = element_text(face = "bold"),
             plot.subtitle = element_text(size = 14, face='bold') )



    # cross-principal plane plots
    cpp_dat_one <- dat %>%
      dplyr::filter(psi == 90)

    cpp_dat <- rbind(cpp_dat_one,
                     cpp_dat_one %>%
                       mutate(tto = -tto))

    cpp <- ggplot(cpp_dat, aes(x = tto, y = BRDF)) +
      geom_line(linetype = "dashed") +
      xlab("VZA") +
      ylab("BRDF") +
      theme_bw() +
      ggtitle("Cross-principal plane") +
      scale_x_continuous(limits = c(-90, 90), breaks = seq(-90, 90, by = 45), expand = c(0, 0)) +
      theme( text = element_text(size = 14, face='bold'),
             axis.title = element_text(size = 16, face = "bold"),
             axis.text = element_text(face = "bold"),
             plot.title = element_text(face = "bold"),
             plot.subtitle = element_text(size = 14, face='bold') )

      gridExtra::grid.arrange(brdf, pp, cpp, layout_matrix = rbind(c(1, 1), c(1, 1), c(1, 1), c(2, 3)))



    showNotification("BRDF plot done successfully.", type = "message")
  })

  # PROSPECT-PRO model table data
  output$prospect_pro_table_brf <- renderTable({
    data.frame(
      Input = c("N", "Cab", "Car", "Anth", "Cbrown", "EWT", "LMA", "alpha", "Prot", "CBC"),
      Description = c(
        "Leaf mesophyll structure index", "Chlorophyll content", "Carotenoid content",
        "Anthocyanins content", "Brown pigments", "Leaf water content",
        "Leaf matter content", "Alpha parameter", "Leaf proteins", "Carbon-based constituent"
      ),
      Units = c("-", "μg cm-2", "μg cm-2", "μg cm-2", "-", "g cm-2", "g cm-2", "-", "g cm-2", "g cm-2"),
      Min = c(1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00),
      Max = c(4.00, 100.00, 40.00, 7.00, 1.00, 0.05, 0.05, 60.00, 0.05, 0.05),
      Default = c(2.50, 40.00, 10.00, 0.00, 0.00, 0.01, 0.00, 40.00, 0.01, 0.01)
    )
  }, rownames = FALSE)

  # fourSAILH model table data
  output$foursailh_table_brf <- renderTable({
    data.frame(
      Input = c("LAI", "LIDFa", "LIDFb", "Hotspot", "tts", "tto", "psi", "psoil"),
      Description = c(
        "Leaf Area Index", "Leaf inclination distribution function a",
        "Leaf inclination distribution function b", "Hot Spot parameter",
        "Sun zenith angle", "Observer zenith angle", "Azimuth Sun / Observer angle",
        "Soil factor"
      ),
      Units = c("m2 m-2", "deg", "deg", "-", "deg", "deg", "deg", "-"),
      Min = c(0.0001, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00),
      Max = c(10.00, 90.00, 0.00, 1.00, 90.00, 90.00, 180.00, 1.00),
      Default = c(4.00, 30.00, 0.00, 0.50, 0.00, 30.00, 0.00, 0.10)
    )
  }, rownames = FALSE)

  # Leaf inclination distribution function (LIDF) table data
  output$lidf_table <- renderTable({
    data.frame(
      Distribution = c("Erectophile", "Planophile", "Spherical", "Plagiophile", "Extremophile","Uniform"),
      Description = c(
        "Leaf angle distribution predominantly upright",
        "Leaf angle distribution predominantly horizontal",
        "Leaf angle distribution is random (spherical)",
        "Leaf angle distribution predominantly inclined",
        "Leaf angle distribution is at extreme angles",
        "Leaf angle distribution is at uniform angles"
      ),
      LIDFa = c(-1, 1, -0.35, 0.00, 0.00, 0.00),
      LIDFb = c(0.00, 0.00, -0.15, -1.00, 1.00, 0.00)
    )
  }, rownames = FALSE)

  # Render data table for the selected dataset
  output$dataset_LUT <- DT::renderDataTable({
    dataset <- reflectance_data_prosail()[[4]]
      DT::datatable(
        dataset,
        options = list(
          scrollX = TRUE,
          columnDefs = list(list(visible = TRUE, targets = 0:4)),
          pageLength = 5  # Set to 10 for viewing 10 rows
        ),
        selection = 'none'
      )

  })


}

# Create Shiny object
shinyApp(ui = ui, server = server)

