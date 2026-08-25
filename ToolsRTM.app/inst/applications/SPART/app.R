
rm(list= ls())

# 1. load the main libraries  -----

library(ToolsRTM)
if (!require(SCOPEinR)) {
  stop("The 'SCOPEinR' package is required. Please install it and try again.\n",
       "To install, run:\n", 
       "  install.packages('SCOPEinR')") 
}
# Required packages
required_packages <- c("shiny", "shinythemes", "shinydashboard", "shinyWidgets", "ggplot2", "dplyr","tidyverse", "reshape2")

# Check for missing packages and install them if necessary
missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

# Load the libraries
lapply(required_packages, library, character.only = TRUE)


# Define UI

ui <- navbarPage("Online reflectance simulator",theme = shinytheme("flatly"),

                 
                 tabPanel("SPART model",
                          # Add a sidebar layout
                          sidebarLayout(
                            # Add a sidebar panel
                            sidebarPanel(
                              # Add a little information about ebird data
                              class = "sidebar",
                              style = "height: 90vh; overflow-y: auto;",
                              # Leaf Model selection
                              h3(''),
                              div(id = "upper-panel",
                                  # Canopy Model selection
                                  selectInput("atmos", label = "Atmospheric Model:",
                                              choices = c("SPART"))
                              ),
                              div(id = "middle-panel",
                                  selectInput("leaf_spart", label = "Leaf Model:",
                                              choices = c("PROSPECT-PRO"))
                                  
                              ),
                              
                              div(id = "lower-panel",
                                  # Canopy Model selection
                                  selectInput("canopy_spart", label = "Canopy Model:",
                                              choices = c("fourSAILH"))
                              ),
                              
                              h4('Atmospheric parameters :'),
                              Parameter = c('Pa','aot550', 'uo3','uh2o','alt_m','Pa0'),
                              Description = c('Air pressure', 'AOT at 550 nm', 'Ozone content ','Water vapour','altitude','sea level air pressure'),
                              # Leaf Model parameters (1)
                              condition = "input.atmos == 'SPART'",
                              sliderInput("Pa_spart", "Air pressure", min = 400, max = 1300, value = 1000, step = 10),
                              sliderInput("aot550_spart", "AOT at 550 nm", min = 0, max = 2.5, value = 0.3246, step = 0.1),
                              sliderInput("uo3_spart", "Ozone content", min = 0, max = 5, value = 0.3480, step = 0.1),
                              sliderInput("uh2o_spart", "Water vapour", min = 0, max = 5, value = 1.4116, step = 0.1),
                              
                              #  sliderInput("alt_m_spart", "altitude", min = 0, max = 600, value = 100, step = 10),
                              #sliderInput("Pa0_spart", "sea level air pressure", min = 999, max = 1020, value = 1013.25, step = 0.1)
                              
                              h4('canopy parameters :'),
                              # Leaf Model parameters
                              conditionalPanel(
                                condition = "input.canopy_spart == 'fourSAILH'",
                                sliderInput("LAI_spart",  HTML("leaf area index (m<sup>2</sup>/m<sup>-2</sup>)"), min = 0, max = 8, value = 4, step = 0.1),
                                sliderInput("LIDFa_spart", "LIDFa (deg)", min = 0.10, max = 90, value = 0.5, step = 0.1),
                                sliderInput("hotspot_spart", "hotspot", min = 0, max = 1, value = 0.5, step = 0.01),
                                sliderInput("tts_spart", "tts (deg)", min = 0, max = 90, value = 0, step = 0.1),
                                sliderInput("tto_spart", "tto (deg)", min = 0, max = 90, value = 30, step = 0.2),
                                sliderInput("psi_spart", "psi (deg)", min = 0, max = 180, value = 45, step = 0.5)
                                
                              ),
                              h4('leaf parameters :'),
                              # Leaf Model parameters (1)
                              
                              conditionalPanel(
                                condition = "input.leaf_spart == 'PROSPECT-PRO'",
                                sliderInput("Cab_spart",HTML("Chlorophyll content (μg cm<sup>-2</sup>)"), min = 0.01, max = 100, value = 50, step = 0.1),
                                sliderInput("Car_spart", HTML("Carotenoids content (μg cm<sup>-2</sup>)"), min = 0, max = 10, value = 2.5, step = 0.1),
                                sliderInput("Anth_spart", HTML("Anthocyanin content (μg cm<sup>-2</sup>)"), min = 0, max = 7, value = 2, step = 0.1),
                                sliderInput("Cbrown_spart", "Cbrown", min = 0, max = 1, value = 0.2, step = 0.1),
                                sliderInput("N_spart", "mesophyll structure parameter", min = 1, max = 4, value = 2.5, step = 0.1),
                                sliderInput("EWT_spart", HTML("Water content (g cm<sup>-2</sup>)"), min = 0.001, max = 0.05, value = 0.01, step = 0.005),
                                sliderInput("Prot_spart",HTML("Proteins (g cm<sup>-2</sup>)"), min = 0.0001, max = 0.03, value = 0.012, step = 0.005),
                                sliderInput("CBC_spart", HTML("Carbon-based constituent (g cm<sup>-2</sup>)"), min = 0, max = 0.03, value = 0.010, step = 0.005)
                                
                              ),
                              
                              #    tags$img(src = "images/JRC.png", height = "100px", width = "200px"),
                              
                            ),# Close sidebarPanel
                            
                            mainPanel(
                              tabsetPanel(
                                ### Firs Pannel
                                tabPanel("Reflectance",
                                         
                                         h2("Interactive reflectance"),
                                         plotOutput("reflectance_plot_spart"),
                                         
                                         
                                         # strong("LUT for leaf + canopy model"),
                                         
                                         # materialSwitch(inputId = "checkbox_scope", label = "Accumulate reflectance spectra", status = "danger"),
                                         
                                         selectInput("sensor_spart", "Select satellite sensor:",
                                                     choices = c('Sentinel2a','Sentinel2b',
                                                                 'MODIS','Landsat-4','Landsat-5','Landsat-7','Landsat-8')),
                                         
                                         p('Download the reflectance spectrum at selected sensor resolution with main plant traits:'),
                                         downloadButton("downloadData_spart", "Download reflectance"),
                                         br(),
                                         br(),
                                         # Conditional panel to show tables after the plot is rendered
                                         conditionalPanel(
                                           condition = "output.reflectance_plot_spart != null", # Show after plot is rendered
                                           p('The selected reflectance spectra were generated based on plant traits and canopy parameters showed in the following tables: '),
                                           
                                           p("Plant traits values used for the selected leaf model:"),
                                           # tableOutput("lut_spart_output_leaf"),
                                           DT::dataTableOutput("lut_spart_output_leaf"),
                                           
                                           p("Structural, viewing angles and geometric values used for the selected canopy model:"),
                                           #tableOutput("lut_spart_output_canopy"),
                                           DT::dataTableOutput("lut_spart_output_canopy"),
                                           strong("Atmospheric parameters:"),
                                           DT::dataTableOutput("lut_spart_output_atmo")
                                           #tableOutput("lut_spart_output_atmo")
                                         )
                                         
                                         
                                         #tableOutput("lut_table.sim")
                                ),         # end Tab Interactive panel
                                
                                tabPanel("Functions",
                                         # Introduction to the module
                                         h4("Main funtions integrated in the ToolsRTM Package"),
                                         
                                         p(style = "text-align: justify;", "This module leverages essential functions from the ", strong("ToolsRTM"), " package to simulate and analyze vegetation data. It particularly focuses on the ", strong("SPART"), " model, which integrates soil-plant-atmosphere radiative transfer for vegetation monitoring using satellite data."),
                                         
                                         # Section for SPART
                                         h4("1. SPART Model"),
                                         p(style = "text-align: justify;", "The ", code("SPART"), "function allow us to simulate the SPART model stands for Soil-Plant-Atmosphere Radiative Transfer, a computationally efficient model for simulating satellite measurements across the solar spectrum. It uses three sub-models: ", strong("BSM (soil)"), ", ", strong("PROSAIL (vegetation canopy)"), ", and ", strong("SMAC (atmosphere)"), " and couples them using the four-stream theory and adding method."),
                                         p(style = "text-align: justify;", "SPART accurately simulates Top of Atmosphere (TOA) spectral observations, considering all major effects like sun-observer geometries and non-Lambertian surface reflectance. This model is particularly useful for analyzing Sentinel-2 imagery and other satellite data, providing accurate simulations of canopy reflectance and atmosphere interaction."),
                                         #br(),
                                         
                                         # Key parameters explanation
                                         p(strong("Key Parameters:")),
                                         # List of function arguments
                                         tags$ul(
                                           tags$li(strong("inputLUT:"), " A Look-Up Table (LUT) containing the distribution of biophysical parameters for modeling."),
                                           tags$li(strong("optipar:"), " Optical parameters and reflection indices from radiative transfer models (RTMs)."),
                                           tags$li(strong("CanopyModel:"), " The canopy model used; default is 'fourSAIL'."),
                                           tags$li(strong("LeafModel:"), " The leaf model to be applied; 'PROSPECT-PRO' is commonly used."),
                                           tags$li(strong("sensor.i:"), " The sensor specification for data alignment; available options include 'Sentinel2A.MSI', 'Sentinel2B.MSI', 'Sentinel3A.OLCI', 'Sentinel3B.OLCI', 'LANDSAT4.TM', 'LANDSAT5.TM', 'LANDSAT7.ETM', 'LANDSAT8.OLI', 'TerraAqua.MODIS'."),
                                           tags$li(strong("df.irradiance:"), " A dataframe of direct and diffuse irradiance values for clear conditions; if NULL, default values will be used."),
                                           tags$li(strong("get.plots:"), " A boolean indicating whether to generate plots; set to TRUE to visualize the results.")
                                         ),
                                         p(strong("Supported sensors:")),
                                         # List of supported sensors
                                         tags$ul(
                                           tags$li(strong("Sentinel2A.MSI:"), " Sentinel-2A MultiSpectral Instrument."),
                                           tags$li(strong("Sentinel2B.MSI:"), " Sentinel-2B MultiSpectral Instrument."),
                                           tags$li(strong("Sentinel3A.OLCI:"), " Sentinel-3A Ocean and Land Colour Instrument."),
                                           tags$li(strong("Sentinel3B.OLCI:"), " Sentinel-3B Ocean and Land Colour Instrument."),
                                           tags$li(strong("LANDSAT4.TM:"), " Landsat 4 Thematic Mapper."),
                                           tags$li(strong("LANDSAT5.TM:"), " Landsat 5 Thematic Mapper."),
                                           tags$li(strong("LANDSAT7.ETM:"), " Landsat 7 Enhanced Thematic Mapper."),
                                           tags$li(strong("LANDSAT8.OLI:"), " Landsat 8 Operational Land Imager."),
                                           tags$li(strong("TerraAqua.MODIS:"), " MODIS onboard Terra and Aqua satellites.")
                                         ),
                                         p(strong("Supported Optical parameters :")),
                                         tags$ul(
                                           tags$li(strong("optipar2017.ProspectD:"), " For the PROSPECT-PRO model."),
                                           tags$li(strong("optipar:"), " For the PROSPECT-D model."),
                                           tags$li(strong("optipar2021.Pro.CX:"), " For the FLUSPECT-B-Cx model."),
                                           tags$li(strong("optipar2020.prospectD.BSM2019:"), " For the FLUSPECT-B model with BSM 2019 parameters.")
                                         ),
                                         
                                         p(strong("Usage:")),
                                         # Example usage of the SPART function
                                         p(style = "text-align: justify;", "Example of using the ", code("SPART"), " function:"),
                                         code("data.spart <- SPART(inputLUT = LUT[1,], optipar = optipar2017.ProspectD"),
                                         br(),
                                         code("CanopyModel = 'fourSAIL',LeafModel = 'PROSPECT-PRO',"),
                                         br(),
                                         code("df.irradiance = NULL,sensor.i = 'Sentinel2A.MSI'  "),
                                         br(),
                                         
                                         # Reference to SPART model documentation
                                         p("More information is available at "),
                                         a("SPART GitHub", href = "https://github.com/peiqiyang/SPART"),
                                         br(),
                                         # Section for SPART
                                         h4("Citation"),
                                         p(style = "text-align: justify;", em("Yang et al. (2020). SPART: Soil-Plant-Atmosphere Radiative Transfer Model for remote sensing applications.")),
                                         
                                         
                                         p(style = "text-align: justify;",HTML('If you use this model with <b>ToolsRTM</b> package, please cite the following references:')),
                                         
                                         p(style = "text-align: justify;",HTML('Camino et al., (2024). RT-Simulator: An Online Platform to Simulate Canopy Reflectance from Biochemical and Structural Plant Properties Using Radiative Transfer Models,
                                          <i>IGARSS 2024 - 2024 IEEE International Geoscience and Remote Sensing Symposium</i>, Athens, Greece, 2024, pp. 2811-2814,
                                          <a href="https://ieeexplore.ieee.org/document/10642442" target="_blank">doi: 10.1109/IGARSS53475.2024.10642442</a>.')),
                                         
                                         p(style = "text-align: justify;",HTML('Arano et al., (2024). Enhancing Chlorophyll Content Estimation With Sentinel-2 Imagery: A Fusion of Deep Learning and Biophysical Models,
                                          <i>IGARSS 2024 - 2024 IEEE International Geoscience and Remote Sensing Symposium</i>, Athens, Greece, 2024, pp. 4486-4489,
                                          <a href="https://ieeexplore.ieee.org/document/10641613" target="_blank">doi: 10.1109/IGARSS53475.2024.10641613</a>.')),
                                         
                                         p(style = "text-align: justify;",'Camino et al., (in prep). Integrating physiological plant traits with
                                             Sentinel-2 imagery for monitoring gross primary production and detecting forest disturbances. '),
                                         
                                         
                                ), #end TabPanel
                                
                                tabPanel("About SPART", "",
                                         wellPanel( style = "background: white",
                                                    h4('ToolsRTM package'),
                                                    p(style = "text-align: justify;", HTML('<b>ToolsRTM</b> package is also designed for running the soil-plant-atmosphere radiative transfer (SPART) model (Yang et. al, 2020) for satellite measurements in the solar spectrum.')),
                                                    
                                                    p(style = "text-align: justify;", 'The SPART model uses three computationally efficient RTMs for soil (BSM), vegetation canopies (PROSAIL) and atmosphere (SMAC), respectively. The sub-models are coupled by using the four-stream theory and the adding method.'),
                                                    p(style = "text-align: justify;", 'The resulting `Soil-Plant-Atmosphere Radiative Transfer model simulates directional TOA spectral observations, with all major effects included, such as sun-observer geometries and non-Lambertia reflectance of the land surface.'),
                                                    
                                                    
                                                    p("More information is available at "),
                                                    a("SPART Github", href = "https://github.com/peiqiyang/SPART"),
                                                    br(),
                                                    h4('SPART model: '),
                                                    
                                                    tableOutput("spart_table"),
                                                    br(),
                                                    
                                                    tags$figure(
                                                      tags$img(src = "images/rtm_spart.png",height = "450px", width = "500px"),
                                                      tags$figcaption("", style = "font-weight: bold;")
                                                      
                                                    ),
                                                    br(),
                                                    strong("Fig 1."),'Simulated reflectance spectrum using the SPART model.' ,
                                                    
                                                    br(),
                                                    h4("Install ToolsRTM"),
                                                    
                                                    p('ToolsRTM is avalaible on gitlab, so you can install using the R console:'),
                                                    
                                                    code('install.packages("toolsrtm-main.tar.gz",repos = NULL,type = "source")'),
                                                    
                                                    
                                                    # Add information on cranes and prompt user to use slider
                                                    
                                                    h4("GitLab repositories"),
                                                    
                                                    a("ToolsRTM package", href = "https://gitlab.com/caminoccg/toolsrtm"),
                                                    br(),
                                                    
                                                    h4("Citation"),
                                                    
                                                    p(style = "text-align: justify;",HTML('If you use <b>ToolsRTM</b> or <b>SCOPEinR</b> packages, please cite the following references:')),
                                                    
                                                    p(style = "text-align: justify;",HTML('Camino et al., (2024). RT-Simulator: An Online Platform to Simulate Canopy Reflectance from Biochemical and Structural Plant Properties Using Radiative Transfer Models,
                                          <i>IGARSS 2024 - 2024 IEEE International Geoscience and Remote Sensing Symposium</i>, Athens, Greece, 2024, pp. 2811-2814,
                                          <a href="https://ieeexplore.ieee.org/document/10642442" target="_blank">doi: 10.1109/IGARSS53475.2024.10642442</a>.')),
                                                    
                                                    p(style = "text-align: justify;",HTML('Arano et al., (2024). Enhancing Chlorophyll Content Estimation With Sentinel-2 Imagery: A Fusion of Deep Learning and Biophysical Models,
                                          <i>IGARSS 2024 - 2024 IEEE International Geoscience and Remote Sensing Symposium</i>, Athens, Greece, 2024, pp. 4486-4489,
                                          <a href="https://ieeexplore.ieee.org/document/10641613" target="_blank">doi: 10.1109/IGARSS53475.2024.10641613</a>.')),
                                                    
                                                    p(style = "text-align: justify;",'Camino et al., (in prep). Integrating physiological plant traits with
                                             Sentinel-2 imagery for monitoring gross primary production and detecting forest disturbances. '),
                                                    
                                                    
                                                    h4("Citation for SPART model"),
                                                    
                                                    p('Yang, P., van der Tol, C., Yin, T., & Verhoef, W. (2020). The SPART model: A soil-plant-atmosphere radiative transfer model for satellite
                                          measurements in the solar spectrum. Remote Sensing of Environment, 247, 111870.'),
                                                    
                                                    p('For the details of the radiative transfer modelling'),
                                                    
                                                    p('Yang, P., Verhoef, W., & van der Tol, C. (2017). The mSCOPE model: A simple adaptation to the SCOPE model to describe reflectance,
                                           fluorescence and photosynthesis of vertically heterogeneous canopies.Remote sensing of environment, 201, 1-11.'),
                                                    #    div(tags$img(src = "images/JRC.png", height = "100px", width = "300px",align="left")
                                                    #       , style="text-align: left;"),
                                                    
                                                    
                                                    
                                                    br(),
                                                    br(),
                                                    br(),
                                                    br(),
                                                    br(),
                                                    
                                         ) # end well panel
                                         
                                         #textInput("Install_package", "", value = "...."),
                                         
                                ), # end Tab Reference panel
                                
                              ) # close tabsetPanel
                              
                            ) # Close mainPanel
                          ) # Close sideBarLayout
                 ), # end tab panel
                 
)
# Define server logic required to draw a histogram ----
server <- function(input, output,session) {
  
  ########### 2) SPART -------------------------------
  
  ### Define table for SPART model ---------------------------------------------------
  
  
  output$spart_table <- renderTable({
    data <- data.frame(
      Parameter = c('Pa','aot550', 'uo3','uh2o','alt_m','Pa0'),
      Description = c('Air pressure', 'AOT at 550 nm', 'Ozone content ','Water vapour','altitude','sea level air pressure'),
      Units = c('hPa','-','cm','g cm-2','m','hPa'),
      Min = c(400,0,0,0,-1000,900),
      Max = c(1300,2.5,5,5,600,1100),
      Default = c(-999,0.3246,0.3480,1.4116,0.0,1013.25)
    )
    colnames(data) <- c("Input", "Description", "Units", "Min", "Max","default")
    return(data)
  })
  
  
  # Reactive expression for the leaf and canopy parameters -------------------------------
  params_spart <- reactive({
    leaf_params <- switch(input$leaf_spart,
                          
                          "PROSPECT-PRO" = c(N=input$N_spart,Cab=input$Cab_spart, Car=input$Car_spart,
                                             Anth=input$Anth_spart, Cbrown=input$Cbrown_spart,
                                             EWT=input$EWT_spart, LMA=0, alpha=40, #input$alpha,
                                             Prot=input$Prot_spart, CBC=input$CBC_spart))
    ## parameters for fourSAIL
    
    canopy_params <- switch(input$canopy_spart,
                            "fourSAILH" = c(LAI=input$LAI_spart, TypeLidf=2, LIDFa=input$LIDFa_spart, LIDFb=0, hspot=input$hotspot_spart,
                                            tts=input$tts_spart, tto=input$tto_spart, psi=input$psi_spart,
                                            LAIu=0.100,sd=640, h=20, cd=4.50,
                                            psoil=0.5))
    
    
    atmos_params <- switch(input$atmos,
                           "SPART" = c(Pa=input$Pa_spart, aot550=input$aot550_spart, uo3=input$uo3_spart,
                                       uh2o=input$uh2o_spart, alt_m=0, Pa0=0.00100133,
                                       skyl=0.1))
    
    
    return(list(leaf_params, canopy_params,atmos_params))
  })
  
  
  # Reactive expression for the LUT data ---------------------------------------------------
  lut_data.leaf.spart <- reactive({
    
    # Get the selected leaf and canopy parameters from the reactive `params` expression
    leaf_params <- unlist(params_spart()[[1]])
    leaf_params_names<-names(params_spart()[[1]])
    
    # Call the LUT function using the selected parameters
    lut_leaf <- data.frame(rbind(leaf_params)) #, canopy_params))
    print(lut_leaf)
    return(lut_leaf)
  })
  
  # Reactive expression for the LUT data ---------------------------------------------------
  lut_data.canopy.spart <- reactive({
    # Get the selected leaf and canopy parameters from the reactive `params` expression
    canopy_params <- unlist(params_spart()[[2]])
    # Call the LUT function using the selected parameters
    lut_canopy <- data.frame(rbind(canopy_params)) #, canopy_params))
    print(lut_canopy)
    return(lut_canopy)
  })
  
  # Reactive expression for the LUT data ---------------------------------------------------
  lut_data.atmos.spart <- reactive({
    # Get the selected leaf and canopy parameters from the reactive `params` expression
    atmos_params <- unlist(params_spart()[[3]])
    # Call the LUT function using the selected parameters
    lut_atmos <- data.frame(rbind(atmos_params)) #, canopy_params))
    print(lut_atmos)
    return(lut_atmos)
  })
  
  # Reactive expression for the LUT data ---------------------------------------------------
  lut_data.sim_spart <- reactive({
    # Get the selected leaf and canopy parameters from the reactive `params` expression
    leaf_params <- unlist(params_spart()[[1]])
    
    canopy_params <- unlist(params_spart()[[2]])
    
    atmos_params <- unlist(params_spart()[[3]])
    
    # Call the LUT function using the selected parameters
    lut_leaf <- data.frame(rbind(leaf_params)) #, canopy_params))
    lut_canopy <- data.frame(rbind(canopy_params)) #, canopy_params))
    lut_atmos <- data.frame(rbind(atmos_params)) #, canopy_params))
    lut_to_sim <-cbind(lut_leaf,lut_canopy,lut_atmos)
    
    # Call the Simulations
    data <- dataSpec_PDB
    Rsoil.dry  <- data[,11]  # rsoil1 = dry soil
    Rsoil.wet <- data[,12]  # rsoil2 = wet soil
    psoil = lut_to_sim[1,'psoil']
    #psoil	 <-  1    # soil factor (psoil=0: wet soil / psoil=1: dry soil)
    rsoil<- c(psoil*Rsoil.dry+(1-psoil)*Rsoil.wet)
    #print(lut_to_sim)
    #reflectance_values <- -foursail(inputLUT=lut.to_sim[1,],rsoil=rsoil,LeafModel = 'PRO')
    return(list(lut_to_sim,rsoil,lut_leaf,lut_canopy,lut_atmos))
    
    
  })
  
  # Reactive expression for the reflectance data ---------------------------------------------------
  reflectance_data_spart <- reactive({
    
    LUT_ <-as.data.frame(lut_data.sim_spart()[[1]])
    
    rsoil_ <- lut_data.sim_spart()[[2]]
    LUT_$TypeLidf = 1
    print(LUT_)
    if (input$sensor_spart == 'Sentinel2a') {
      
      sensor.i = Sentinel2A.MSI
      
      
    } else if (input$sensor_spart == 'Sentinel2b') {
      
      sensor.i = Sentinel2B.MSI
      
      
    } else if (input$sensor_spart == 'Landsat-4') {
      
      sensor.i = LANDSAT4.TM
      
    } else if (input$sensor_spart == 'Landsat-5') {
      
      sensor.i = LANDSAT5.TM
      
      
    } else if (input$sensor_spart == 'Landsat-7') {
      
      sensor.i = LANDSAT7.ETM
      
      
    } else if (input$sensor_spart == 'Landsat-8') {
      
      sensor.i = LANDSAT8.OLI
      
      
    } else if (input$sensor_spart == 'MODIS') {
      sensor.i = TerraAqua.MODIS
      
    }
    
    data.spart<- SPART(inputLUT = LUT_[1,],optipar=optipar2021.Pro.CX,
                       CanopyModel = 'fourSAIL',
                       LeafModel='PROSPECT-PRO',
                       df.irradiance = NULL,
                       sensor.i = sensor.i)
    
    rfl.toa <- data.spart$output$rfl.toa
    rfl.toc <- data.spart$output$rfl.toc
    rfl.toc.brdf <- data.spart$output$rfl.toc.BRDF
    
    reflectance_df <- data.frame(wavelength = data.spart$output$wave, rfl.toa = rfl.toa,
                                 rfl.toc=rfl.toc,rfl.toc.brdf=rfl.toc.brdf)
    showNotification("Satellite selection with applied bandwidths completed successfully.", type = "message")
    
    return(reflectance_df)
    
  })
  
  
  # Render the reflectance plot ---------------------------------------------------
  output$reflectance_plot_spart <- renderPlot({
    
    
    to_plot<-reflectance_data_spart()
    print(to_plot)
    ggplot(data = to_plot, aes(x = wavelength)) +
      labs(y = "Reflectance", x = "") +
      geom_point(aes(y = rfl.toa, color = "TOA rfl."), size = 1) +
      geom_line(aes(y = rfl.toa, color = "TOA rfl."),size = 1.0) +
      
      geom_point(aes(y = rfl.toc, color = "TOC rfl. (SMAC)"), size = 1.0) +
      geom_line(aes(y = rfl.toc, color = "TOC rfl. (SMAC)"), size = 1.0) +
      
      theme_bw() +
      guides(color = guide_legend(title = ""), linetype = guide_legend(title = ""), shape = guide_legend(title = "")) +
      
      theme(legend.position = "top",
            legend.key.size = unit(4, "lines"),
            text = element_text(size = 16, face='bold'),  # Increase the text size
            legend.text = element_text(face = "bold", size = 14), # Increase legend text size and make it bold
            axis.title.x = element_text(face = "bold", size = 14), # Increase x axis label size and make it bold
            axis.title.y = element_text(face = "bold", size = 14)) # Increase y axis label size and make it bold
    
    
    
  })
  
  # Create a reactive data frame for saving data
  save_data_spart <- reactive({
    to_save<-reflectance_data_spart()
  })
  # Add the download handler for saving data ------
  output$downloadData_spart <- downloadHandler(
    filename = function() {
      paste("SPART_sim_with_", input$leaf_spart,'_', input$canopy_spart,'_', input$sensor_spart, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(save_data_spart(), file, row.names = FALSE)
    }
  )
  
  
  #### Render the LUT leaf table (SPART) ----------------------------------
  
  output$lut_spart_output_leaf <- DT::renderDataTable({
    DT::datatable(lut_data.sim_spart()[[3]],
                  options = list(
                    dom = 't',    # Removes search panel, pagination, etc.
                    pageLength = 1 # Shows only 1 row
                  ), rownames = FALSE  # Removes row names
    )
  })
  
  
  #### Render the LUT canopy table (SPART) ----------------------------------
  
  output$lut_spart_output_canopy <- DT::renderDataTable({
    DT::datatable(lut_data.sim_spart()[[4]],
                  options = list(
                    dom = 't',    # Removes search panel, pagination, etc.
                    pageLength = 1 # Shows only 1 row
                  ), rownames = FALSE  # Removes row names
    )
  })
  
  #### Render the LUT atmospheric table (SPART) ----------------------------------
  
  output$lut_spart_output_atmo <- DT::renderDataTable({
    DT::datatable(lut_data.sim_spart()[[5]],
                  options = list(
                    dom = 't',    # Removes search panel, pagination, etc.
                    pageLength = 1 # Shows only 1 row
                  ), rownames = FALSE  # Removes row names
    )
  })
  
  
  

}

# Create Shiny object
shinyApp(ui = ui, server = server)

