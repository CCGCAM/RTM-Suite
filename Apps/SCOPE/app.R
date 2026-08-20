
rm(list= ls())

# 1. load the main libraries  -----

library(ToolsRTM)
if (!require(SCOPEinR)) {
  stop("The 'SCOPEinR' package is required. Please install it and try again.\n",
       "To install, run:\n", 
       "  install.packages('SCOPEinR')") 
}
required_packages <- c("shiny", "shinythemes", "shinybusy", 'shinyWidgets',"ggplot2", "dplyr",
                       "doParallel",'foreach','parallel','DT')

# Check for missing packages and install them if necessary
missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}
# Load the libraries
lapply(required_packages, library, character.only = TRUE)

# Define UI

ui <- navbarPage("Online reflectance simulator",theme = shinytheme("flatly"),
                 
                 tabPanel("SCOPE model",
                    
                          # Add a sidebar layout
                          sidebarLayout(
                            # Add a sidebar panel
                            sidebarPanel(
                              # Add a little information about ebird data
                              class = "sidebar",
                              style = "height: 90vh; overflow-y: auto;",
                              # Leaf Model selection
                              h3('Select RT models:'),
                              div(id = "upper-panel",
                                  selectInput("Leaf_scope", label = "Leaf Model:",
                                              choices = c("FLUSPECT-B", "FLUSPECT-B-Cx"))
                                  
                              ),
                              div(id = "lower-panel",
                                  # Canopy Model selection
                                  selectInput("Canopy_scope", label = "Canopy Model:",
                                              choices = c("fourSAILH"))
                              ),
                              
                              h4('leaf parameters :'),
                              # Leaf Model parameters (1)
                              
                              conditionalPanel(
                                condition = "input.Leaf_scope == 'FLUSPECT-B'",
                                sliderInput("fqe_fd_scope", "Fluorescence quantum efficiency (fqe)", min = 0.01, max = 0.05, value = 0.02, step = 0.005),
                                sliderInput("Cab_fd_scope",  HTML("Chlorophyll content (μg cm<sup>-2</sup>)"), min = 0.01, max = 100, value = 50, step = 0.1),
                                sliderInput("Car_fd_scope",  HTML("Carotenoid content (μg cm<sup>-2</sup>)"), min = 0.01, max = 40, value = 20, step = 0.1),
                                sliderInput("Cs_fd_scope", "Leaf Senescence", min = 0, max = 1, value = 0.01),
                                sliderInput("Cx_fd_scope", "Violaxanthin - Zeaxanthin transition status", min = 0, max = 1, value = 0.1),
                                
                                sliderInput("EWT_fd_scope",  HTML("Water content (g cm<sup>-2</sup>)"), min = 0.001, max = 0.05, value = 0.01, step = 0.005),
                                sliderInput("LMA_fd_scope", HTML("Dry matter content (g cm<sup>-2</sup>)"), min = 0, max = 0.05, value = 0.005),
                                sliderInput("N_fd_scope", "mesophyll structure parameter", min = 1, max = 4, value = 2.5, step = 0.1)
                                
                              ),
                              
                              # Leaf Model parameters (4)
                              conditionalPanel(
                                condition = "input.Leaf_scope == 'FLUSPECT-B-Cx'",
                                sliderInput("fqe_fp_scope", "Fluorescence quantum efficiency (fqe)", min = 0.01, max = 0.05, value = 0.02, step = 0.005),
                                sliderInput("Cab_fp_scope",  HTML("Chlorophyll content (μg cm<sup>-2</sup>)"), min = 0.01, max = 100, value = 50, step = 0.1),
                                sliderInput("Car_fp_scope",  HTML("Carotenoid content (μg cm<sup>-2</sup>)"), min = 0.01, max = 40, value = 20, step = 0.1),
                                sliderInput("Anth_fp_scope", HTML("Anthocyanin content (μg cm<sup>-2</sup>)"), min = 0.01, max = 7, value = 2, step = 0.1),
                                sliderInput("Cs_fp_scope", "Leaf Senescence", min = 0, max = 1, value = 0.01),
                                sliderInput("Cx_fp_scope", "Violaxanthin - Zeaxanthin transition status", min = 0, max = 1, value = 0.1),
                                sliderInput("EWT_fp_scope", HTML("Water content (g cm<sup>-2</sup>)"), min = 0, max = 0.05, value = 0.01, step = 0.005),
                                sliderInput("LMA_fp_scope",  HTML("Dry matter content (g cm<sup>-2</sup>)"), min = 0, max = 0.05, value = 0.005),
                                sliderInput("N_fp_scope", "mesophyll structure parameter", min = 1, max = 4, value = 2.5, step = 0.1)
                              ),
                              
                              
                              h4('leaf biochemical parameters :'),
                              # Leaf Model parameters (1)
                              
                              sliderInput("Vcmax_scope", HTML("Vcmax (μmol m<sup>-2</sup>s<sup>-1</sup>)"), min = 0, max = 175, value = 80, step = 0.5),
                              sliderInput("BaLBerrySlope", "Ball-Berry Slope", min = 0, max = 8, value = 8, step = 0.1),
                              sliderInput("BaLBerry0", "BaLBerry 0", min = 0, max = 0.2, value = 0.1, step = 0.001),
                              
                              h4('Meteorological parameters :'),
                              # Leaf Model parameters (1)
                              
                              sliderInput("Rn_scope",  HTML("broadband incoming shortwave radiation (Watt m<sup>-2</sup>)"), min = 400, max = 1000, value = 600, step = 10),
                              sliderInput("Rli_scope", HTML("broadband incoming longwave radiation (Watt m<sup>-2</sup>)"), min = 150, max = 450, value = 300, step = 10),
                              sliderInput("Ta_scope", HTML("Air temperature (deg)"), min = 20, max = 40, value = 25, step = 0.1),
                              
                              
                              h4('canopy parameters :'),
                              # Leaf Model parameters
                              conditionalPanel(
                                condition = "input.Canopy_scope == 'fourSAILH'",
                                sliderInput("LAI_scope",  HTML("leaf area index (m<sup>2</sup>/m<sup>-2</sup>)"), min = 0, max = 8, value = 4, step = 0.1),
                                sliderInput("LIDFa_scope", "LIDFa", min = -1, max = 1, value = 0.5, step = 0.1),
                                sliderInput("LIDFb_scope", "LIDFb", min = -1, max = 1, value = 0.5, step = 0.1),
                                sliderInput("hotspot_scope", "hotspot", min = 0, max = 1, value = 0.5, step = 0.01),
                                sliderInput("tts_scope", "tts (deg)", min = 0, max = 90, value = 0, step = 0.1),
                                sliderInput("tto_scope", "tto (deg)", min = 0, max = 90, value = 30, step = 0.2)),
                              
                              #    tags$img(src = "images/JRC.png", height = "100px", width = "200px"),
                              
                            ),# Close sidebarPanel
                            
                            mainPanel(
                              tabsetPanel(
                                ### Firs Pannel
                                tabPanel("Reflectance",
                                         
                                         h2("Interactive reflectance"),
                                         
                                         plotOutput("reflectance_plot_scope"),
                                         # strong("LUT for leaf + canopy model"),
                                         
                                         materialSwitch(inputId = "checkbox_scope", label = "Accumulate reflectance spectra", status = "danger"),
                                         
                                         selectInput("sensor_scope", "Select satellite sensor:",
                                                     choices = c('SCOPE','Sentinel2a','Sentinel2b',
                                                                 'EnMAP','Hyperion','MODIS',
                                                                 'Landsat4','Landsat5','Landsat7','Landsat8', #'ALI',
                                                                 'Quickbird','RapidEye','WorldView2-4','WorldView2-8')),
                                         
                                         p('Download the reflectance spectrum at selected sensor resolution with main plant traits:'),
                                         downloadButton("downloadData_scope", "Download reflectance"),
                                         br(),
                                         br(),
                                         p('The selected reflectance spectra were generated based on plant traits and canopy parameters showed in the following tables: '),
                                         br(),
                                         br(),
                                         #uiOutput("plotReady_tableReady_scope"),
                                         strong("Plant traits values used for the selected leaf model:"),
                                         br(),
                                         br(),
                                         #verbatimTextOutput("lut_scope_output_leaf"),
                                         # Output for the DataTable
                                         DT::dataTableOutput("lut_table_scope_leaf"),
                                         strong("Plant traits and geometric parameters values used for the selected canopy model:"),
                                         br(),
                                         br(),
                                         DT::dataTableOutput("lut_table_scope_canopy"),
                                         #verbatimTextOutput("lut_scope_output_canopy"),
                                         br(),
                                         br(),
                                         DT::dataTableOutput("lut_table_scope_bioleaf"),
                                         
                                         #tableOutput("lut_table.sim")
                                ),         # end Tab Interactive panel
                                
                                ### Second Pannel
                                tabPanel("Fluorescence",
                                         h2("Chlorohyl fluorescence emission"),
                                         
                                         plotOutput("sif_plot_scope"),
                                         
                                         materialSwitch(inputId = "checkbox_sif", label = "Accumulate SIF spectra", status = "danger"),
                                         
                                         selectInput("sensor_sif", "Select satellite sensor:",
                                                     choices = c('SCOPE','Sentinel2a','Sentinel2b',
                                                                 'EnMAP','Hyperion','MODIS',
                                                                 'Landsat4','Landsat5','Landsat7','Landsat8', #'ALI',
                                                                 'Quickbird','RapidEye','WorldView2-4','WorldView2-8')),
                                         p('Download the reflectance spectrum at selected sensor resolution with main plant traits:'),
                                         downloadButton("downloadData_sif", "Download fluorescence"),
                                         br(),
                                         br(),
                                         p('The selected fluorescence spectra were generated based on plant traits and canopy parameters showed in the Reflectance tab: '),
                                         br(),
                                         br(),
                                         
                                ),# end Tab Inputs panel
                                # end Tab Inputs panel
                                tabPanel("Radiance",
                                         h2("Radiance profile"),
                                         plotOutput("Lo_plot_scope"),
                                         
                                         materialSwitch(inputId = "checkbox_LO", label = "Accumulate radiance spectra", status = "danger"),
                                         
                                         selectInput("sensor_LO", "Select satellite sensor:",
                                                     choices = c('SCOPE','Sentinel2a','Sentinel2b',
                                                                 'EnMAP','Hyperion','MODIS',
                                                                 'Landsat4','Landsat5','Landsat7','Landsat8', #'ALI',
                                                                 'Quickbird','RapidEye','WorldView2-4','WorldView2-8')),
                                         p('Download the radiance spectrum at selected sensor resolution with main plant traits:'),
                                         downloadButton("downloadData_LO", "Download radiance"),
                                         br(),
                                         br(),
                                         p('The selected radiance spectra were generated based on plant traits and canopy parameters showed in the Reflectance tab: '),
                                         br(),
                                         br(),
                                         # Add a verbatimTextOutput for printing LUT.SCOPE
                                         #verbatimTextOutput("lut_scope_output"),
                                         br(),
                                         #verbatimTextOutput("model_scope_outputs"),
                                         br(),
                                         #      verbatimTextOutput("lut_scope_output_leaf"),
                                         
                                         br(),
                                         #     verbatimTextOutput("lut_scope_output_param"),
                                         
                                         br(),
                                         
                                         
                                ),# end Tab Reference panel
                                tabPanel("Functions",
                                         # Introduction to the module
                                         h4("Main funtions integrated in the SCOPEinR Package"),
                                         
                                         p(style = "text-align: justify;", "This module utilizes key functions integrated into the ", strong("SCOPEinR"), " package to process and analyze vegetation data derived from Sentinel-2 satellite imagery. These functions are critical for estimating spectral indices, training machine learning models, and predicting plant traits such as gross primary production (GPP) and other physiological traits."),
                                         
                                         # getLUTfromRanges
                                         h4("1. getLUT.SCOPE"),
                                         # Section for getMLmodel.withRetrain
                                         p(style = "text-align: justify;", "The ", code("getLUT.SCOPE"), " function generates a Look-Up Tables (LUTs) adaptaed to SCOPE model. It requires a table of inputs specifying the ranges of variables and allows users
                      to define the number of LUT rows and set a seed for randomization. This function is essential for creating LUTs used in various analyses related to vegetation and canopy properties.
                      Below are the key arguments for this function:"),
                                         ## Key parameters explanation
                                         p(strong("Key Parameters:")),
                                         # List of function arguments
                                         tags$ul(
                                           tags$li(strong("inputs:"), " A table with specific ranges and types of distributions for the variables."),
                                           tags$li(strong("nLUT:"), " The number of rows to be generated for the LUT."),
                                           tags$li(strong("setseed:"), " A seed number to control the random process. By default, it is set to 123, but you can specify a different number to ensure reproducibility.")
                                         ),
                                         p(strong("Usage:")),
                                         # Example usage of the function
                                         p(style = "text-align: justify;", "Example of using this function:"),
                                         code("LUTs <- getLUT.SCOPE(inputs = SCOPEinR::inputsSCOPE, nLUT = 100, setseed = 123);"),
                                         br(),
                                         # Section for getSCOPE.parallel
                                         h4("2. get.SCOPE.parallel"),
                                         p(style = "text-align: justify;", "The ", code("get.SCOPE.parallel"), " function runs SCOPE simulations in parallel, allowing for efficient processing of large datasets. This function requires a Look-Up Table (LUT) containing various properties needed for the SCOPE model, and it offers options to customize the output, parallel processing, and the models used. Below are the key arguments for this function:"),
                                         
                                         # Key parameters explanation
                                         p(strong("Key Parameters:")),
                                         # List of function arguments
                                         tags$ul(
                                           tags$li(strong("LUT:"), " Leaf, biochemistry, viewing angles, meteo, and canopy properties needed for running the SCOPE model."),
                                           tags$li(strong("options.SCOPE:"), " Optical leaf properties and total irradiance. If a specific value is provided, it will replace the usual Modtran output."),
                                           tags$li(strong("path.out:"), " Folder for saving the SCOPE outputs."),
                                           tags$li(strong("parallel:"), " Logical, indicating whether to use parallel processing. Default is TRUE."),
                                           tags$li(strong("canopy.model:"), " Selection of canopy model options available: 'fourSAIL', 'INFORM'. Default is 'fourSAIL'."),
                                           tags$li(strong("leaf.model:"), " Selection of leaf model options available: 'fluspect-CX', 'fluspect-B', 'PROSPECT', 'Liberty'. Default is 'fluspect-CX'."),
                                           tags$li(strong("get.outputs:"), " Specifies which variables to retrieve; 'ALL' retrieves all variables, while 'Main' retrieves only the main variables."),
                                           tags$li(strong("get.plots:"), " Logical, indicating whether to plot intermediate results. Default is TRUE."),
                                           tags$li(strong("get.csv:"), " Logical, indicating whether to save outputs. Default is TRUE."),
                                           tags$li(strong("n.cores:"), " Integer, indicating the number of cores to use. Default is 2 if this parameter is null or missing.")
                                         ),
                                         
                                         p(strong("Usage:")),
                                         p(strong("Usage:")),
                                         p(style = "text-align: justify;", "Example of using this function:"),
                                         code("db.sims <- SCOPEinR::get.SCOPE.parallel(LUT = my_LUT, options.SCOPE = my_options,"),
                                         br(),
                                         code(" path.out = 'outputs/', parallel = TRUE, canopy.model = 'fourSAIL', "),
                                         br(),
                                         code("leaf.model = 'fluspect-CX',"),
                                         br(),
                                         code("get.outputs = 'Main', get.plots = TRUE, get.csv = TRUE, n.cores = 8)"),
                                         
                                         h4("Citation "),
                                         p(style = "text-align: justify;",HTML('If you use these functions with <b>SCOPEinR</b> package, please cite the following references:')),
                                         
                                         p(style = "text-align: justify;",HTML('Camino et al., (2024). RT-Simulator: An Online Platform to Simulate Canopy Reflectance from Biochemical and Structural Plant Properties Using Radiative Transfer Models,
                                          <i>IGARSS 2024 - 2024 IEEE International Geoscience and Remote Sensing Symposium</i>, Athens, Greece, 2024, pp. 2811-2814,
                                          <a href="https://ieeexplore.ieee.org/document/10642442" target="_blank">doi: 10.1109/IGARSS53475.2024.10642442</a>.')),
                                         
                                         p(style = "text-align: justify;",HTML('Arano et al., (2024). Enhancing Chlorophyll Content Estimation With Sentinel-2 Imagery: A Fusion of Deep Learning and Biophysical Models,
                                          <i>IGARSS 2024 - 2024 IEEE International Geoscience and Remote Sensing Symposium</i>, Athens, Greece, 2024, pp. 4486-4489,
                                          <a href="https://ieeexplore.ieee.org/document/10641613" target="_blank">doi: 10.1109/IGARSS53475.2024.10641613</a>.')),
                                         
                                         p(style = "text-align: justify;",'Camino et al., (in prep). Integrating physiological plant traits with
                                             Sentinel-2 imagery for monitoring gross primary production and detecting forest disturbances. '),
                                         
                                ),
                                tabPanel("About SCOPEinR", "",
                                         wellPanel( style = "background: white",
                                                    h4('SCOPEinR package'),
                                                    p(style = "text-align: justify;",HTML('<b>SCOPEinR</b> package is designed for running the Soil Canopy Observation, Photochemistry and Energy fluxes (SCOPE, Van der Tol at al., 2009, Yang et al., 2020) radiative transfer model.')),
                                                    
                                                    p(style = "text-align: justify;",'This package simulates reflectance and chorophyll fluorescence emission using the SCOPE model. To make inter-comparison with other main radiative transfer (RT) models is recommeded to install the ToolsRTM package. This ToolsRTM package uses several functions for simulating canopy reflectance at several spectral resolution (1nm, hyper spectral and Sentinel-2).'),
                                                    
                                                    "Manual is available at ",
                                                    tags$a(href="https://carlos-camino.shinyapps.io/0-toolsrtm-simulator/_w_e851c6b2/Notebooks/R/SCOPEinR/SCOPEinR.html",
                                                           "ReadTheDocs"),
                                                    ".",
                                                    br(),
                                                    
                                                    tags$figure(
                                                      tags$img(src = "images/rtm_scope.png",height = "450px", width = "500px"),
                                                      tags$figcaption("", style = "font-weight: bold;")
                                                      
                                                    ),
                                                    br(),
                                                    strong("Fig 1."),'Simulated reflectance spectrum using the SCOPE model.' ,
                                                    #tags$img(src = "rtm_sims.png", height = "450px", width = "500px"),
                                                    #htmltools::img(src = "rtm_sims.png"),
                                                    
                                                    br(),
                                                    br(),
                                                    tags$figure(
                                                      tags$img(src = "images/rtm_scope_sif.png",height = "450px", width = "500px"),
                                                      tags$figcaption("", style = "font-weight: bold;")
                                                      
                                                    ),
                                                    br(),
                                                    strong("Fig 1."),'Simulated chorohyll fluorescence emission using theSCOPE model.' ,
                                                    
                                                    br(),
                                                    h4("Install SCOPEinR"),
                                                    p('SCOPEinR is avalaible on gitlab, so you can install using the R console:'),
                                                    
                                                    code('install.packages("scopeinr-main.tar.gz",repos = NULL,type = "source")'),
                                                    
                                                    
                                                    # Add information on cranes and prompt user to use slider
                                                    
                                                    h4("GitLab repositories"),
                                                    
                                                    a("SCOPEinR package", href = "https://gitlab.com/caminoccg/scopeinr"),
                                                    br(),
                                                    h3("Citation"),
                                                    
                                                    p(HTML('If you use <b>ToolsRTM</b> or <b>SCOPEinR</b> packages, please cite the following references:')),
                                                    
                                                    p(style = "text-align: justify;",HTML('Camino et al., (2024). RT-Simulator: An Online Platform to Simulate Canopy Reflectance from Biochemical and Structural Plant Properties Using Radiative Transfer Models,
                                           <i>IGARSS 2024 - 2024 IEEE International Geoscience and Remote Sensing Symposium</i>, Athens, Greece, 2024, pp. 2811-2814,
                                          <a href="https://ieeexplore.ieee.org/document/10642442" target="_blank">doi: 10.1109/IGARSS53475.2024.10642442</a>.')),
                                                    
                                                    p(style = "text-align: justify;",HTML('Arano et al., (2024). Enhancing Chlorophyll Content Estimation With Sentinel-2 Imagery: A Fusion of Deep Learning and Biophysical Models,
                                          <i>IGARSS 2024 - 2024 IEEE International Geoscience and Remote Sensing Symposium</i>, Athens, Greece, 2024, pp. 4486-4489,
                                          <a href="https://ieeexplore.ieee.org/document/10641613" target="_blank">doi: 10.1109/IGARSS53475.2024.10641613</a>.')),
                                                    
                                                    p(style = "text-align: justify;",'Camino et al., (in prep). Integrating physiological plant traits with
                                             Sentinel-2 imagery for monitoring gross primary production and detecting forest disturbances. '),
                                                    
                                                    h4("Citation for SCOPE model"),
                                                    
                                                    p(style = "text-align: justify;",'Yang, P., Prikaziuk, E., Verhoef, W., and Van der Tol, C. 2021 "SCOPE 2.0: a model to simulate vegetated land surface
                                         fluxes and satellite signals" Geoscientific Model Development, 14, 4697–4712, https://doi.org/10.5194/gmd-14-4697-2021'),
                                                    
                                                    p('Van der Tol, C., W. Verhoef, J Timmermans, A Verhoef, and Z Su. 2009. “An Integrated Model of Soil-Canopy Spectral Radiances,
                                         Photosynthesis, Fluorescence, Temperature and Energy Balance.” Biogeosciences 6 (12): 3109–29, https://doi.org/10.5194/bg-6-3109-2009'),
                                                    
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
  
  ########### 3) SCOPE -------------------------------
  
  ## Reactive expression (SCOPE) for the leaf and canopy parameters -------------------------------
  params_SCOPE <- reactive({
    leaf_params_SCOPE <- switch(input$Leaf_scope,
                                
                                
                                "FLUSPECT-B" = c(N=input$N_fd_scope,Cab = input$Cab_fd_scope, Car = input$Car_fd_scope,
                                                 EWT=input$EWT_fd_scope,LMA=input$LMA_fd_scope,
                                                 Cx=input$Cx_fd_scope, Cs = input$Cs_fd_scope,
                                                 fqe = input$fqe_fd_scope),
                                
                                "FLUSPECT-B-Cx" = c(N=input$N_fp_scope,Cab = input$Cab_fp_scope, Car = input$Car_fp_scope,
                                                    EWT=input$EWT_fp_scope,LMA=input$LMA_fp_scope,
                                                    Cx=input$Cx_fp_scope, Cs = input$Cs_fp_scope,
                                                    fqe = input$fqe_fp_scope))
    
    
    ## // parameters for fourSAIL
    
    
    canopy_params_SCOPE <- switch(input$Canopy_scope,
                                  "fourSAILH" = c(LAI=input$LAI_scope, TypeLidf=1, LIDFa=input$LIDFa_scope,
                                                  LIDFb=input$LIDFb_scope, hspot=input$hotspot_scope,
                                                  tts=input$tts_scope, tto=input$tto_scope,
                                                  psi=0.5))
    
    
    ## // dditional parameters
    
    Vcmax_scope <- input$Vcmax_scope
    BaLBerrySlope <- input$BaLBerrySlope
    BaLBerry0 <- input$BaLBerry0
    
    Rn_scope <- input$Rn_scope
    Rli_scope <- input$Rli_scope
    Ta_scope <- input$Ta_scope
    
    additional_params_SCOPE <- c(Vcmax25 = Vcmax_scope, BallBerrySlope = BaLBerrySlope,
                                 BallBerry0 = BaLBerry0,
                                 Rin = Rn_scope,
                                 Rli = Rli_scope,
                                 Ta =Ta_scope)
    
    return(list(leaf_params_SCOPE, canopy_params_SCOPE,additional_params_SCOPE))
  })
  
  
  # Reactive expression for the LUT leaf data (SCOPE)---------------------------------------------------
  lut_scope.leaf <- reactive({
    
    # Get the selected leaf and canopy parameters from the reactive `params` expression
    leaf_params <- unlist(params_SCOPE()[[1]])
    leaf_params_names<-names(params_SCOPE()[[1]])
    
    bioleaf_params <- unlist(params_SCOPE()[[3]])
    bioleaf_params_names<-names(params_SCOPE()[[3]])
    # Call the LUT function using the selected parameters
    
    # Call the LUT function using the selected parameters
    lut_leaf <- c(leaf_params, bioleaf_params)
    #colnames(lut_leaf) <- c(leaf_params_names, names(params_SCOPE()[[3]]))
    
    # lut_leaf <- data.frame(rbind(leaf_params, bioleaf_params_names)) #, canopy_params))
    print(lut_leaf)
    return(lut_leaf)
  })
  
  
  
  # Reactive expression for LUT canopy data ---------------------------------------------------
  lut_scope.canopy <- reactive({
    # Get the selected leaf and canopy parameters from the reactive `params` expression
    canopy_params <- unlist(params_SCOPE()[[2]])
    # Call the LUT function using the selected parameters
    lut_canopy <- canopy_params
    return(lut_canopy)
  })
  
  
  
  # Render the printed output in verbatimTextOutput
  #output$lut_scope_output_param <- renderPrint({
  # params_SCOPE()
  
  #})
  
  #### Render the LUT SCOPE param  ----------------------------------
  
  output$lut_scope_output_param <- DT::renderDataTable({
    DT::datatable(params_SCOPE(),
                  options = list(
                    dom = 't',    # Removes search panel, pagination, etc.
                    pageLength = 1 # Shows only 1 row
                  ), rownames = FALSE  # Removes row names
    )
  })
  
  
  # Render the printed output in verbatimTextOutput
  output$lut_scope_output_canopy <- renderPrint({
    lut_scope.canopy()
  })
  
  # Render the printed output in verbatimTextOutput
  output$lut_scope_output_leaf <- renderPrint({
    lut_scope.leaf()
  })
  
  #### Render the LUT SCOPE canopy  ----------------------------------
  
  output$lut_scope_output_canopy <- DT::renderDataTable({
    DT::datatable(lut_scope.canopy(),
                  options = list(
                    dom = 't',    # Removes search panel, pagination, etc.
                    pageLength = 1 # Shows only 1 row
                  ), rownames = FALSE  # Removes row names
    )
  })
  
  
  
  # Render the printed output in verbatimTextOutput
  output$lut_scope_output.leaf <- renderPrint({
    lut_scope.leaf()
    print(class(lut_scope.leaf()))
    
  })
  
  # Render the LUT SCOPE leaf parameters ----------------------------------
  
  output$lut_table_scope_leaf <- DT::renderDataTable({
    DT::datatable(data.frame(t(unlist(params_SCOPE()[[1]]))),
                  options = list(
                    dom = 't',    # Removes search panel, pagination, etc.
                    pageLength = 1 # Shows only 1 row
                  ), rownames = FALSE  # Removes row names
    )
  })
  
  output$lut_table_scope.leaf <- renderTable({
    lut_scope.leaf()
  })
  
  
  
  
  #### Render the LUT SCOPE for canopy parameters  ----------------------------------
  
  output$lut_table_scope_canopy <- DT::renderDataTable({
    DT::datatable(data.frame(t(unlist(params_SCOPE()[[2]]))),
                  options = list(
                    dom = 't',    # Removes search panel, pagination, etc.
                    pageLength = 1 # Shows only 1 row
                  ), rownames = FALSE  # Removes row names
    )
  })
  
  #### Render the LUT SCOPE for bioleaf parameters  ----------------------------------
  
  output$lut_table_scope_bioleaf <- DT::renderDataTable({
    DT::datatable(data.frame(t(unlist(params_SCOPE()[[3]]))),
                  options = list(
                    dom = 't',    # Removes search panel, pagination, etc.
                    pageLength = 1 # Shows only 1 row
                  ), rownames = FALSE  # Removes row names
    )
  })
  
  
  # Reactive expression for the LUT data ---------------------------------------------------
  lut_scope.sim <- reactive({
    # Get the selected leaf and canopy parameters from the reactive `params` expression
    leaf_params <- unlist(params_SCOPE()[[1]])
    bioleaf_params <- unlist(params_SCOPE()[[3]])
    
    canopy_params <- unlist(params_SCOPE()[[2]])
    
    # Call the LUT function using the selected parameters
    lut_leaf <- data.frame(t(c(leaf_params,bioleaf_params)) )#, canopy_params))
    
    lut_canopy <- data.frame(t(canopy_params)) #, canopy_params))
    
    lut.to_sim <-data.frame(cbind(lut_leaf,lut_canopy))
    
    
    # Call the Simulations
    LUT.SCOPE <- SCOPE.LUT.default
    
    # Assign values from lut.to_sim to LUT.SCOPE
    for (input_name in names(lut.to_sim)) {
      LUT.SCOPE[[input_name]] <- lut.to_sim[[input_name]]
    }
    
    
    
    # Print the LUT.SCOPE dataframe
    print(LUT.SCOPE)
    return(list(LUT.SCOPE))
    
    
  })
  
  
  # Render the printed output in verbatimTextOutput
  output$lut_scope_output <- renderPrint({
    lut_scope.sim()
    
    
  })
  
  
  
  
  
  # Reactive expression for the reflectance data ---------------------------------------------------
  get_scope <- reactive({
    
    LUT_ <- lut_scope.sim()[[1]]
    
    showNotification("Processing the SCOPE model using SCOPEinR package. Please wait...", type = "message")
    showNotification("SCOPE requires additional time due to the model's complexity.", type = "warning")
    
    db.sim <- get.SCOPE(LUT=LUT_,options.SCOPE = data.opts,path.out = 'www/outs/',
                        optipar=optipar2021.Pro.CX,
                        leaf.model='fluspect-CX',canopy.model='fourSAIL',
                        get.outputs = 'Main', get.plots = F)
    showNotification("Simulation completed successfully using the SCOPEinR package.", type = "message")
    
    
    reflectance_app <- db.sim[[1]]$data.rad$reflapp
    reflectance_values <- db.sim[[1]]$data.rad$refl
    reflectance_rsd <- db.sim[[1]]$data.rad$rsd
    
    reflectance_rdo <- db.sim[[1]]$data.rad$rdo
    
    wave.rfl <- db.sim[[1]]$data.spectral$wlS
    
    LoF_sunlit <- db.sim[[1]]$data.rad$LoF_sunlit
    LoF_shaded <- db.sim[[1]]$data.rad$LoF_shaded
    LoF_soil<- db.sim[[1]]$data.rad$LoF_soil
    LoF_all<- db.sim[[1]]$data.rad$Femliave_
    
    LoF_ <- db.sim[[1]]$data.rad$LoF_
    wave.LoF <- db.sim[[1]]$data.spectral$wlF
    
    Lotot_ <- db.sim[[1]]$data.rad$Lotot_
    Lototf_ <- db.sim[[1]]$data.rad$Lototf_
    wave.Lot <- db.sim[[1]]$data.spectral$wlS
    # Combine the reflectance values with wavelength values
    reflectance_df <- data.frame(wavelength = wave.rfl,
                                 reflapp = reflectance_app,
                                 reflectance = reflectance_values,rdo =reflectance_rdo,
                                 rsd = reflectance_rsd)
    
    
    LoF_df <- data.frame(wavelength = wave.LoF,
                         LoF_sunlit = LoF_sunlit, LoF_soil = LoF_soil,
                         LoF_shaded = LoF_shaded, LoF_all = LoF_all,
                         LoF_ = LoF_)
    Lotot_df <- data.frame(wavelength = wave.Lot, Lotot_ = Lotot_,Lototf_ = Lototf_)
    
    
    #  return(list(reflectance_df = reflectance_df, LoF_df= LoF_df, Lotot_df= Lotot_df))
    return(list(reflectance_df = reflectance_df, LoF_df= LoF_df, Lotot_df= Lotot_df))
    
    
  })
  
  
  # Render the printed output in verbatimTextOutput
  output$model_scope_outputs <- renderPrint({
    get_scope()
    
    
  })
  
  
  
  # Render the reflectance plot ---------------------------------------------------
  output$reflectance_plot_scope <- renderPlot({
    
    # Show the modal window
    show_modal_spinner()
    
    df.to_plot<-get_scope()$reflectance_df
    
    # Plot the variables using ggplot
    p1 <- ggplot(df.to_plot, aes(x = wavelength)) +
      geom_rect(aes(xmin = 725, xmax = 800, ymin = 0.3, ymax = 0.7), color = 'grey', linetype = 1, alpha = 0.1) +
      geom_line(aes(y = reflectance, color = 'reflectance'), linewidth = 1) +
      geom_line(aes(y = reflapp, color = 'reflect. Apparent'), linewidth = 1) +
      geom_line(aes(y = rdo, color = 'rdo'), linewidth = 1) +
      geom_line(aes(y = rsd, color = 'rsd'), linewidth = 1) +
      labs(x = "Wavelength (nm)", y = "Reflectance") +
      #    scale_color_manual(values = c('black','orange', 'red', 'navyblue'),
      #                      labels = c('reflectance','reflect.App','rdo', 'rsd')) +
      theme_bw() +
      xlim(400, 2400) + ylim(0,0.7) +
      theme(legend.position = 'top',
            #       legend.box.background = element_rect(color = "black",linewidth=1),
            plot.title = element_text(hjust = 0.5, size=14,face="bold"),
            axis.title = element_text(face="bold", size=14),
            legend.text = element_text(face="bold", size=10),
            axis.text.y=element_text(hjust = 0.5, size=12,face="bold"),
            axis.text.x=element_text(hjust = 0.5, size=12,face="bold"),
            legend.title=element_blank())
    
    # Zoomed-in plot
    p2 <- ggplot(df.to_plot, aes(x = wavelength)) +
      geom_line(aes(y = reflectance, color = 'reflectance'), linewidth = 1) +
      geom_line(aes(y = reflapp, color = 'reflect. Apparent'), linewidth = 1) +
      geom_line(aes(y = rdo, color = 'rdo'), linewidth = 1) +
      geom_line(aes(y = rsd, color = 'rsd'), linewidth = 1) +
      labs(x = "", y = "reflectance") +
      theme_bw() +
      xlim(725, 800) +
      ylim(0.2, 0.7) +
      theme(
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
        axis.title = element_text(face = "bold", size = 10),
        axis.text.y = element_text(hjust = 0.5, size = 8, face = "bold"),
        axis.text.x = element_text(hjust = 0.5, size = 8, face = "bold"),
        legend.title = element_blank()
      )
    
    # Remove modal window when done
    showNotification("SCOPE executed successfully.", type = "message")
    remove_modal_spinner()
    # Combine the plots
    p1 +
      annotation_custom(ggplotGrob(p2), xmin = 410, xmax = 700, ymin = 0.3, ymax = 0.7) +
      geom_rect(aes(xmin = 400, xmax = 710, ymin = 0.3, ymax = 0.7), color = 'black', linetype = 'dashed', alpha = 0)
    
    # geom_path(
    #  aes(x, y, group = grp),
    # data = data.frame(x = c(800, 700, 900, 800), y = c(0, 0.4, 0, 0.7), grp = c(1, 1, 2, 2)),
    #  linetype = 'dashed')
    
    
  })
  
  
  
  # Create a reactive expression for the accumulated data
  accumulated_data <- reactive({
    if (input$checkbox_sif) {
      # Accumulate the data
      accumulate_spectra(get_scope()$LoF_df)
    } else {
      get_scope()$LoF_df
    }
  })
  
  
  # Accumulating spectra -------------------
  accumulate_spectra <- function(df) {
    # //Melt the data frame to long format
    df.melted <- reshape2::melt(df, id.vars = 'wavelength', measure.vars = "LoF_", variable.name = 'Variable', value.name = "Spectrum")
    
    # Repeat the spectra columns based on the number of repeats you want
    num_repeats <- 5  # Specify the number of repeats you want
    df.melted <- df.melted %>%
      group_by(Variable) %>%
      mutate(Spectrum = rep(Spectrum, num_repeats))
    
    # Return the accumulated data frame
    df.melted
  }
  
  
  
  
  # Render the SIF plot ---------------------------------------------------
  output$sif_plot_scope<- renderPlot({
    
    
    
    if (input$checkbox_sif == F) {
      
      
      df.to_plot <- accumulated_data()
      ggplot(df.to_plot, aes(x = wavelength)) +
        geom_line(aes(y = LoF_sunlit, color = 'LoF sunlit'), linewidth = 1) +
        # geom_line(aes(y = LoF_shaded, color = 'LoF shaded'), linewidth = 1) +
        #  geom_line(aes(y = LoF_soil, color = 'LoF soil'), linewidth = 1) +
        #  geom_line(aes(y = LoF_all, color = 'LoF all leaves'), linewidth = 1) +
        geom_line(aes(y = LoF_, color = 'LoF_'), linewidth = 1) +
        
        labs(x = "Wavelength (nm)", y = "fluorescence emission") +
        #   scale_color_manual(values = c('orange', 'red','black'),
        #                     labels = c('LoF_sunlit','LoF_shaded','LoF_total')) +
        theme_bw() + xlim(640,800) +  theme(legend.position = 'top',
                                            #  legend.box.background = element_rect(color = "black",linewidth=1),
                                            plot.title = element_text(hjust = 0.5, size=14,face="bold"),
                                            axis.title = element_text(face="bold", size=14),
                                            legend.text = element_text(face="bold", size=10),
                                            axis.text.y=element_text(hjust = 0.5, size=12,face="bold"),
                                            axis.text.x=element_text(hjust = 0.5, size=12,face="bold"),
                                            legend.title=element_blank()) +
        # legend.title = element_text(face = "bold", size = 14)) +
        guides(color = guide_legend(title.position = "top", title.hjust = 0.5))
    } else {
      
      df.to_plot <- accumulated_data()
      # Plot the spectra for each variable
      ggplot(df.to_plot, aes(x = wavelength, y = Spectrum, color = )) +
        geom_line(linewidth = 1) +
        labs(x = "Wavelength (nm)", y = "fluorescence emission") +
        #    scale_color_manual(values = color_palette1(num_inputs1)) +
        #   scale_color_manual(values = color_palette2(num_inputs2)) +
        theme_bw()
    }
    
  })
  
  
  ## Render the Lototf_ plot ---------------------------------------------------
  output$Lo_plot_scope<- renderPlot({
    
    df.to_plot<-get_scope()$Lotot_df
    
    p1 <-ggplot(df.to_plot, aes(x = wavelength)) +
      geom_rect(aes(xmin = 725, xmax = 800, ymin = 20, ymax = 150), color = 'grey', linetype = 1, alpha = 0.1) +
      
      geom_line(aes(y = Lototf_, color = 'Lototf_'), linewidth = 1) +
      geom_line(aes(y = Lotot_, color = 'Lotot_'), linewidth = 1) +
      
      labs(x = "Wavelength (nm)", y = "Radiance excluding/adding fluorescence") +
      theme_bw()  + xlim(400,900) +
      #    scale_color_manual(values = c('orange', 'black'),
      #                      labels = c('Lototf_','Lotot_')) +
      theme(legend.position = 'top',
            #  legend.box.background = element_rect(color = "black",linewidth=1),
            plot.title = element_text(hjust = 0.5, size=14,face="bold"),
            axis.title = element_text(face="bold", size=14),
            legend.text = element_text(face="bold", size=10),
            axis.text.y=element_text(hjust = 0.5, size=12,face="bold"),
            axis.text.x=element_text(hjust = 0.5, size=12,face="bold"),
            legend.title=element_blank()) +
      # legend.title = element_text(face = "bold", size = 14)) +
      guides(color = guide_legend(title.position = "top", title.hjust = 0.5))
    
    
    # Zoomed-in plot
    p2 <- ggplot(df.to_plot, aes(x = wavelength)) +
      geom_line(aes(y = Lototf_, color = 'Lototf_'), linewidth = 1) +
      geom_line(aes(y = Lotot_, color = 'Lotot_'), linewidth = 1) +
      
      labs(x = "", y = "Radiance") +
      theme_bw() +
      xlim(725, 800) +
      ylim(0, 150) +
      theme(
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
        axis.title = element_text(face = "bold", size = 10),
        axis.text.y = element_text(hjust = 0.5, size = 8, face = "bold"),
        axis.text.x = element_text(hjust = 0.5, size = 8, face = "bold"),
        legend.title = element_blank())
    
    p1 +
      annotation_custom(ggplotGrob(p2), xmin = 410, xmax = 600, ymin = 50, ymax = 150) +
      geom_rect(aes(xmin = 400, xmax = 610, ymin = 50, ymax = 150), color = 'black', linetype = 'dashed', alpha = 0)
    
    
    
  })
  
  ########## end SCOPE module ######################################################
  
}

# Create Shiny object
shinyApp(ui = ui, server = server)

