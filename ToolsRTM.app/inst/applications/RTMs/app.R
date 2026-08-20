
rm(list= ls())

# 1. load the main libraries  -----

library(ToolsRTM)
required_packages <- c("shiny", "shinythemes", 'shinyWidgets',"ggplot2", "dplyr",'DT')

# Check for missing packages and install them if necessary
missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}
# Load the libraries
lapply(required_packages, library, character.only = TRUE)

# Define UI

ui <- navbarPage("Online reflectance simulator",theme = shinytheme("flatly"),
                 
                 tabPanel(title = "Interactive ToolsRTM",
                       
                          # Add a sidebar layout
                          sidebarLayout(
                            # Add a sidebar panel
                            sidebarPanel(
                              # Add a little information about ebird data
                              class = "sidebar",
                              style = "height: 90vh; overflow-y: auto;",
                              # Leaf Model selection
                              h3(''),
                              #tags$img(src= "images/JRC.png", height = 80, width = 200),
                              div(id = "upper-panel",
                                  selectInput("leaf_model", label = "Leaf Model:",
                                              choices = c("PROSPECT-PRO","PROSPECT-D", "Liberty", "FLUSPECT-B", "FLUSPECT-B-Cx"))
                                  
                              ),
                              div(id = "lower-panel",
                                  # Canopy Model selection
                                  selectInput("canopy_model", label = "Canopy Model:",
                                              choices = c("fourSAILH", "INFORM")) # fourSAILH2
                              ),
                              h4('leaf parameters :'),
                              # Leaf Model parameters (1)
                              # Checkbox to show/hide leaf parameters
                              #checkboxInput("show_leaf_params", "Show leaf parameters", value = FALSE),
                              
                              conditionalPanel(
                                condition = "input.leaf_model == 'PROSPECT-PRO'",
                                sliderInput("Cab", HTML("Chlorophyll content (μg cm<sup>-2</sup>)"), min = 0, max = 100, value = 45, step = 0.5),
                                sliderInput("Car", HTML("Carotenoid content (μg cm<sup>-2</sup>)"), min = 0, max = 20, value = 12, step = 0.5),
                                sliderInput("Anth", HTML("Anthocyanin content (μg cm<sup>-2</sup>)"), min = 0, max = 7, value = 0, step = 0.5),
                                sliderInput("Cbrown", "Cbrown", min = 0, max = 1, value = 0.2, step = 0.05),
                                sliderInput("N", "mesophyll structure parameter", min = 1, max = 4, value = 2.45, step = 0.1),
                                sliderInput("EWT", HTML("Water content (g cm<sup>-2</sup>)"), min = 0.0001, max = 0.05, value = 0.01, step = 0.005),
                                #    sliderInput("LMA", "LMA (g cm-2)", min = 0, max = 0.05, value = 0.005),
                                #  sliderInput("alpha", "alpha", min = 0, max = 60, value = 40, step = 0.1),
                                sliderInput("Prot", HTML("Proteins (g cm<sup>-2</sup>)"), min = 0.0001, max = 0.03, value = 0.012, step = 0.005),
                                sliderInput("CBC",  HTML("Carbon-based constituent (g cm<sup>-2</sup>)"), min = 0, max = 0.03, value = 0.010, step = 0.005)
                                
                              ),
                              
                              # Leaf Model parameters (2)
                              conditionalPanel(
                                condition = "input.leaf_model == 'PROSPECT-D'",
                                
                                sliderInput("Cab_d",  HTML("Chlorophyll content (μg cm<sup>-2</sup>)"), min = 0, max = 100, value = 20, step = 0.1),
                                sliderInput("Car_d", HTML("Carotenoid content (μg cm<sup>-2</sup>)"), min = 0, max = 20, value = 2.5, step = 0.1),
                                sliderInput("Anth_d",  HTML("Anthocyanin content (μg cm<sup>-2</sup>)"), min = 0, max = 7, value = 2, step = 0.1),
                                sliderInput("Cbrown_d", "Cbrown", min = 0, max = 1, value = 0.2, step = 0.1),
                                sliderInput("N_d", "mesophyll structure parameter", min = 1, max = 4, value = 2.5, step = 0.1),
                                sliderInput("EWT_d", HTML("Water content (g cm<sup>-2</sup>)"), min = 0, max = 0.05, value = 0.009, step = 0.002),
                                sliderInput("LMA_d", HTML("dry matter content (g cm<sup>-2</sup>)"), min = 0, max = 0.05, value = 0.009, step = 0.005)
                                #  sliderInput("alpha_d", "alpha", min = 0, max = 60, value = 20, step = 0.1)
                                
                              ),
                              # Leaf Model parameters (3)
                              conditionalPanel(
                                condition = "input.leaf_model == 'Liberty'",
                                sliderInput("Cab_l", HTML("Chlorophyll content (μg cm<sup>-2</sup>)"), min = 0, max = 60, value = 40, step = 1),
                                sliderInput("EWT_l", HTML("Water content (g cm<sup>-2</sup>)"), min = 0, max = 0.05, value = 0.009, step = 0.001),
                                sliderInput("lign_cell", "Lignin and cellulose content", min = 10, max = 80, value = 40, step = 0.1),
                                sliderInput("Nitrogen", HTML("Nitrogen content (g m<sup>-2</sup>)"), min = 0.3, max = 2, value = 1, step = 0.1),
                                sliderInput("cell_d", HTML("Cell diameter  (m<sup>-6</sup>)"), min = 20, max = 200, value = 58, step = 0.1),
                                sliderInput("inter_c", "Intercellular air space", min = 0.01, max = 0.1, value = 0.045, step = 0.005),
                                sliderInput("baseline_abs", "baseline", min = 0.0004, max = 0.0006, value = 0.0004),
                                sliderInput("leaf_thick", "leaf thickness", min = 1, max = 10, value = 1.6, step = 0.1),
                                sliderInput("albino_abs", "Albino absorption", min = 0, max = 4, value = 2, step = 0.1)
                                
                                
                              ),
                              # Leaf Model parameters (4)
                              conditionalPanel(
                                condition = "input.leaf_model == 'FLUSPECT-B'",
                                
                                sliderInput("fqe_fd", "Fluorescence quantum efficiency (fqe)", min = 0, max = 0.05, value = 0.02, step = 0.005),
                                sliderInput("Cab_fd", HTML("Chlorophyll content (μg cm<sup>-2</sup>)"), min = 0, max = 100, value = 50, step = 0.1),
                                sliderInput("Car_fd", HTML("Carotenoid content (μg cm<sup>-2</sup>)"), min = 0, max = 20, value = 20, step = 0.1),
                                sliderInput("Cs_fd", "Leaf Senescence", min = 0, max = 1, value = 0.01),
                                sliderInput("Cx_fd", "Violaxanthin - Zeaxanthin transition status", min = 0, max = 1, value = 0.1),
                                sliderInput("N_fd", "mesophyll structure parameter", min = 1, max = 4, value = 2.5, step = 0.1),
                                sliderInput("EWT_fd", HTML("Water content (g cm<sup>-2</sup>)"), min = 0, max = 0.05, value = 0.01, step = 0.005),
                                sliderInput("LMA_fd", HTML("dry matter content (g cm<sup>-2</sup>)"), min = 0, max = 0.05, value = 0.01)
                                
                              ),
                              
                              # Leaf Model parameters (4)
                              conditionalPanel(
                                condition = "input.leaf_model == 'FLUSPECT-B-Cx'",
                                sliderInput("fqe_fp", "Fluorescence quantum efficiency (fqe)", min = 0, max = 0.05, value = 0.02, step = 0.005),
                                sliderInput("Cab_fp", HTML("Chlorophyll content (μg cm<sup>-2</sup>)"), min = 0, max = 100, value = 50, step = 0.1),
                                sliderInput("Car_fp", HTML("Carotenoid content (μg cm<sup>-2</sup>)"), min = 0, max = 40, value = 20, step = 0.1),
                                sliderInput("Anth_fp",  HTML("Anthocyanin content (μg cm<sup>-2</sup>)"), min = 0, max = 7, value = 2, step = 0.1),
                                sliderInput("Cs_fp", "Leaf Senescence", min = 0, max = 1, value = 0.01),
                                sliderInput("Cx_fp", "Violaxanthin - Zeaxanthin transition status", min = 0, max = 1, value = 0.1),
                                sliderInput("N_fp", "mesophyll structure parameter", min = 1, max = 4, value = 2.5, step = 0.1),
                                sliderInput("EWT_fp", HTML("Water content (g cm<sup>-2</sup>)"), min = 0, max = 0.05, value = 0.01, step = 0.005),
                                sliderInput("LMA_fp", HTML("dry matter content (g cm<sup>-2</sup>)"), min = 0, max = 0.05, value = 0.005),
                                sliderInput("Prot_fp", HTML("Proteins (g cm<sup>-2</sup>)"), min = 0, max = 0.03, value = 0.012, step = 0.005),
                                sliderInput("CBC_fp", HTML("Carbon-based Constituent  (g cm<sup>-2</sup>)"), min = 0, max = 0.03, value = 0.010, step = 0.005)
                              ),
                              
                              h4('canopy parameters :'),
                              # Leaf Model parameters
                              conditionalPanel(
                                condition = "input.canopy_model == 'fourSAILH'",
                                sliderInput("LAI", HTML("leaf area index (m<sup>2</sup>/m<sup>-2</sup>)"), min = 0, max = 8, value = 4, step = 0.1),
                                sliderInput("LIDFa", "LIDFa (°)", min = 0, max = 90, value = 30, step = 0.1),
                                sliderInput("hotspot", "hotspot", min = 0, max = 1, value = 0.5, step = 0.01),
                                sliderInput("tts", "tts (deg)", min = 0, max = 90, value = 0, step = 0.1),
                                sliderInput("tto", "tto (deg)", min = 0, max = 90, value = 30, step = 0.2),
                                sliderInput("psi", "psi (deg)", min = 0, max = 180, value = 0, step = 0.5),
                                sliderInput("psoil", "soil factor", min = 0, max = 1, value = 0.35,step = 0.1)
                                
                              ),
                              conditionalPanel(
                                condition = "input.canopy_model == 'INFORM'",
                                sliderInput("LAI_", HTML("leaf area index (m<sup>2</sup>/m<sup>-2</sup>)"), min = 0, max = 8, value = 4, step = 0.1),
                                sliderInput("LIDFa_", "LIDFa (deg)", min = 0, max = 90, value = 30, step = 0.1),
                                sliderInput("hotspot_", "hotspot", min = 0, max = 1, value = 0.1, step = 0.1),
                                sliderInput("tts_", "tts (deg)", min = 0, max = 90, value = 0, step = 0.1),
                                sliderInput("tto_", "tto (deg)", min = 0, max = 90, value = 45, step = 0.2),
                                sliderInput("psi_", "psi (deg)", min = 0, max = 180, value = 0, step = 0.5),
                                sliderInput("LAIu_", HTML("understory LAI (m<sup>2</sup>/m<sup>-2</sup>)"), min = 0.05, max = 2, value = 0.5,step = 0.015),
                                sliderInput("cd_", "Crown diameter (m)", min = 0, max = 10, value = 4.5, step = 0.5),
                                sliderInput("sd_",  HTML("Stem density (ha<sup>-1</sup>)"), min = 0, max = 3000, value = 2500, step = 1.5),
                                sliderInput("h_", "tree height (m)", min = 1, max = 40, value = 20, step = 0.1),
                                sliderInput("skyl_", "skyl (fraction-Fixed to 0.1)", min = 0.001, max = 0.4, value = 0.1, step = 0.005),
                                sliderInput("psoil_", "soil factor", min = 0, max = 1, value = 0.1, step = 0.1)
                              ),
                              
                              conditionalPanel(
                                condition = "input.canopy_model == 'fourSAILH2'",
                                sliderInput("p1_s1", HTML("leaf area index (m<sup>2</sup>/m<sup>-2</sup>)"), min = 0, max = 8, value = 4, step = 0.1),
                                sliderInput("p2_s2", "tts (deg)", min = 0, max = 90, value = 0, step = 0.1),
                                sliderInput("p3_s2", "tto (deg)", min = 0, max = 90, value = 45, step = 0.2),
                                sliderInput("p4_s2", "psi (deg)", min = 0, max = 180, value = 0, step = 0.5)
                              ),
                              # tags$img(src = "images/JRC.png", height = "100px", width = "200px"),
                              
                            ),# Close sidebarPanel
                            
                            mainPanel(
                              tabsetPanel(
                                ### Firs Pannel
                                tabPanel("Canopy reflectance",
                                         h2("Interactive reflectance"),
                                         #checkboxInput("checkbox", "Accumulate reflectance spectra", FALSE),
                                         #prettyCheckbox(inputId = "checkbox", label = "Check me!", icon = icon("check")),
                                         plotOutput("reflectance_plot"),
                                         
                                         materialSwitch(inputId = "checkbox", label = "Accumulate reflectance spectra", status = "danger"),
                                         
                                         selectInput("sensor", "Select satellite sensor:",
                                                     choices = c('RTM','Sentinel2a','Sentinel2b',
                                                                 'EnMAP','PRISMA','Hyperion','MODIS',
                                                                 'Landsat4','Landsat5','Landsat7','Landsat8', #'ALI',
                                                                 'Quickbird','RapidEye','WorldView2-4','WorldView2-8')),
                                         checkboxInput("AddRegion", "Adding bandwidths", value = FALSE),
                                         p('Download the reflectance spectrum at selected sensor resolution with main plant traits:'),
                                         downloadButton("downloadData_sim", "Download spectrum"),
                                         br(),
                                         br(),
                                         p('The selected reflectance spectrum was generated based on plant traits and canopy parameters showed in the following tables: '),
                                         
                                         strong("Plant traits values used for the selected leaf model:"),
                                         br(),
                                         br(),
                                         DT::dataTableOutput("lut_table.leaf"),
                                         # tableOutput("lut_table.leaf"),
                                         strong("Plant traits and geometric parameters values used for the selected canopy model:"),
                                         br(),
                                         br(),
                                         DT::dataTableOutput("lut_table.canopy"),
                                         #tableOutput("lut_table.canopy"),
                                         # strong("LUT for leaf + canopy model"),
                                         br(),
                                         br(),
                                         #tableOutput("lut_table.sim")
                                         #DT::dataTableOutput("lut_table.sim"),
                                ),         # end Tab Interactive panel
                                
                                ### Second Pannel
                                tabPanel("Inputs",
                                         h4("PROSPECT model"),
                                         p(style = "text-align: justify;",'Plant leaf reflectance and transmittance are calculated from 400 nm to 2500 nm (1 nm step) with the following parameters:'),
                                         tableOutput("prospect_table"),
                                         br(),
                                         p(style = "text-align: justify;",HTML('<b>PROSPECT-D </b>model uses the following parameters: Cab, Car, Anth, CBrown, alpha, EWT and LMA')),
                                         p(style = "text-align: justify;",HTML('<b>PROSPECT-PRO </b>model uses the following parameters: Cab, Car, Anth, CBrown, alpha, EWT, leaf proteins and CBC')),
                                         p('In PROSPECT-PRO model, LMA = 0, when leaf proteins and CBC is provided'),
                                         
                                         h4("Liberty model"),
                                         
                                         p(style = "text-align: justify;",'Leaf radiative transfer model designed for conifer needles'),
                                         p(style = "text-align: justify;",HTML('<b>Liberty </b>model simulates plant leaf reflectance and transmittance  with the following parameters:')),
                                         
                                         tableOutput("Liberty_table"),
                                         br(),
                                         p(style = "text-align: justify;",'Note that absorption coefficients used in Liberty model are based on the original values from Dawson et al. (1998). Improved specific absorption coefficients are available from later work by Di Vittorio (2009).'),
                                         
                                         h4("FLUSPECT model"),
                                         
                                         p(style = "text-align: justify;",'Leaf radiative transfer model designed for adding fluorescence emision.'),
                                         tableOutput("fluspect_table"),
                                         h4("fourSAILH model"),
                                         
                                         p(HTML('<b>fourSAILH </b>model is based on a version provided by	Wout Verhoef et al. (2007)')),
                                         p("original version downloadable at ",
                                           tags$a(href="http://teledetection.ipgp.jussieu.fr/prosail/", "http://teledetection.ipgp.jussieu.fr/prosail/")),
                                         tableOutput("foursailh_table"),
                                         h4("INFORM model"),
                                         
                                         p(style = "text-align: justify;",HTML('<b>INFORM </b>model simulates the bi-directional reflectance of forest stands between 400 and 2500 nm. ')),
                                         p(style = "text-align: justify;",HTML('This model integrates the  Ground coverage <b>(FLIM)</b> model and computes the average horizontal area of a single tree crown in hectare (k) corrected by the factor (adapt=0.6)')),
                                         tableOutput("inform_table"),
                                         
                                         
                                         
                                ),# end Tab Inputs panel
                                tabPanel("Functions",
                                         # Introduction to the module
                                         h4("Main funtions integrated in the ToolsRTM Package"),
                                         
                                         p(style = "text-align: justify;", "This module utilizes key functions integrated into the ", strong("ToolsRTM"), " package to process and analyze vegetation data derived from Sentinel-2 satellite imagery. These functions are critical for estimating spectral indices, training machine learning models, and predicting plant traits such as gross primary production (GPP) and other physiological traits."),
                                         
                                         # Section for fourSAIL
                                         h4("1. fourSAIL"),
                                         p(style = "text-align: justify;", "The ", code("fourSAIL"), " function models canopy reflectance using various leaf models such as ", strong("PROSPECT-D"), ", ", strong("PROSPECT-PRO"), ", ", strong("Liberty"), ", ", strong("FLUSPECT-B-Cx"), ", and ", strong("FLUSPECT-B"), ". It calculates the Bidirectional Reflectance Function (BRF) based on vegetation parameters, soil reflectance, and canopy structure, providing accurate reflectance simulations at Sentinel-2 resolution."),
                                         # Key parameters explanation
                                         p(strong("Key Parameters:")),
                                         # List of function arguments
                                         tags$ul(
                                           tags$li(strong("inputLUT:"), " A Look-Up Table (LUT) containing the distribution of biophysical parameters used as input parameters in the model."),
                                           tags$li(strong("rsoil:"), " A numeric value representing the soil reflectance, which affects the overall canopy reflectance calculation."),
                                           tags$li(strong("LeafModel:"), " The version of the leaf model used; options include 'PROSPECT-PRO', 'PROSPECT-D', 'Liberty', 'FLUSPECT-B' and 'FLUSPECT-B-Cx'. The default is 'PROSPECT-PRO'."),
                                           tags$li(strong("spectrum.all:"), " A boolean value indicating whether to calculate the full spectrum. Set to TRUE for PROSPECT and Liberty models (400-2500 nm) and FALSE for SPART and Fluspect models (400-2400 nm).")
                                         ),
                                         p(strong("Usage:")),
                                         # Example usage of the function
                                         p(style = "text-align: justify;", "Example of using the ", code("fourSAIL"), " function:"),
                                         code("reflectance <- fourSAIL(inputLUT = myLUTTable, rsoil = 0.2, LeafModel = 'PROSPECT-PRO', spectrum.all = TRUE);"),
                                         br(),
                                         
                                         
                                         # Section for INFORM
                                         h4("2. INFORM"),
                                         p(style = "text-align: justify;", "The ", code("INFORM"), " function models forest canopy reflectance, with a focus on analyzing disturbances such as bark beetle outbreaks. It integrates spectral data from leaf and soil models to assess forest health. This function simulates reflectance across the visible, near-infrared, and shortwave infrared wavelengths for various platforms such as Sentinel-2 and Landsat."),
                                         p(style = "text-align: justify;", em("Reference: Camino et al. (2025). Understanding bark beetle outbreaks using INFORM canopy models in forest ecosystems.")),
                                         
                                         # Key parameters explanation
                                         p(strong("Key Parameters:")),
                                         # List of function arguments
                                         tags$ul(
                                           tags$li(strong("inputLUT:"), " A Look-Up Table (LUT) that specifies the distribution of biophysical parameters used as input parameters in the model."),
                                           tags$li(strong("rsoil:"), " A numeric value representing the soil reflectance, which influences the overall reflectance calculation."),
                                           tags$li(strong("LeafModel:"), " The version of the leaf model to be used; options include 'PROSPECT-PRO','PROSPECT-D', 'Liberty', 'FLUSPECT-B' and 'FLUSPECT-B-Cx'."),
                                         ),
                                         p(strong("Usage:")),
                                         # Example usage of the function
                                         p(style = "text-align: justify;", "Example of using the ", code("INFORM"), " function:"),
                                         code("reflectance <- INFORM(inputLUT = LUT_[1,], rsoil = rsoil_, LeafModel = 'PROSPECT-PRO');"),
                                         br(),
                                         
                                         # Section for Compute_BRF
                                         h4("3. Compute_BRF"),
                                         p(style = "text-align: justify;", "The ", code("Compute_BRF"), " function computes the Bidirectional Reflectance Function (BRF) using shortwave infrared (SWIR) spectral data. It models the interaction of light with the canopy, incorporating both soil and leaf scattering properties. This function provides detailed simulations across the spectrum from visible to shortwave infrared."),
                                         p(style = "text-align: justify;", "This function is particularly beneficial for models that require SWIR wavelengths (up to 2400 nm), as utilized in ", code("FLUSPECT-B"), " and ", code("FLUSPECT-B-Cx"), "."),
                                         
                                         # Key parameters explanation for Compute_BRF
                                         p(strong("Key Parameters:")),
                                         tags$ul(
                                           tags$li(strong("rdot:"), " Hemispherical-directional reflectance factor in the viewing direction, representing how much light is reflected by the canopy in a specified direction."),
                                           tags$li(strong("rsot:"), " Bi-directional reflectance factor, which represents the reflectance of the soil surface from multiple angles."),
                                           tags$li(strong("tts:"), " A numeric value indicating the viewing angle (top of the canopy)."),
                                           tags$li(strong("data.light:"), " A data frame containing the  direct solar radiation and diffuse solar radiation required for the BRF calculation; If this is NULL, we use a default values")
                                         ),
                                         
                                         p(strong("Usage:")),
                                         # Example usage of the Compute_BRF function
                                         p(style = "text-align: justify;", "Example of using the ", code("Compute_BRF"), " function:"),
                                         code("reflectance_values <- foursail(inputLUT = LUT_[1,], rsoil = rsoil_, LeafModel = 'PROSPECT-PRO');"),
                                         br(),
                                         code("rdot <- reflectance_values[[1]];"),
                                         br(),
                                         code("rsot <- reflectance_values[[2]];"),
                                         br(),
                                         code("brf_values <- Compute_BRF(rdot = rdot, rsot = rsot, tts = LUT_[1, 'tts'], data.light = dataSpec_PDB);"),
                                         br(),
                                         
                                         # Section for get.spectral.convolution.rfl
                                         h4("4. get.spectral.convolution.rfl"),
                                         p(style = "text-align: justify;", "The ", code("get.spectral.convolution.rfl"), " function resamples reflectance data to align with specific satellite bands (e.g., Landsat-8, Landsat-5, Sentinel-2). It takes high-resolution reflectance data and convolves it to match the spectral response of various sensors, making the data compatible for different satellite imagery applications."),
                                         p(style = "text-align: justify;", "This function ensures that reflectance values from simulations are accurately aligned with the spectral bands of the satellite platform."),
                                         
                                         # Key parameters explanation for get.spectral.convolution.rfl
                                         p(strong("Key Parameters:")),
                                         tags$ul(
                                           tags$li(strong("df:"), " A data frame containing the high-resolution reflectance data to be resampled."),
                                           tags$li(strong("sensor.i:"), " The sensor specification to which the data will be aligned."),
                                           tags$li(strong("get.plots:"), " A boolean value indicating whether to generate plots of the convolution results (default is FALSE).")
                                         ),
                                         
                                         
                                         # Section for Sensor Specification
                                         h4("Sensor Specification (sensor.i)"),
                                         p(style = "text-align: justify;", "In ToolsRTM, the available sensor specifications for aligning reflectance data include:"),
                                         
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
                                         
                                         p(strong("Usage:")),
                                         # Example usage of the get.spectral.convolution.rfl function
                                         p(style = "text-align: justify;", "Example of using the ", code("get.spectral.convolution.rfl"), " function:"),
                                         code("df_resampled <- get.spectral.convolution.rfl(df = Spec.simula, sensor.i = ToolsRTM::Sentinel2A.MSI, get.plots = FALSE);"),
                                         br(),# Section for get.spectral.convolution.rfl
                                         h4("Citation"),
                                         p(style = "text-align: justify;",HTML('If you use these models with <b>ToolsRTM</b> package, please cite the following references:')),
                                         
                                         p(style = "text-align: justify;",HTML('Camino et al., (2024). RT-Simulator: An Online Platform to Simulate Canopy Reflectance from Biochemical and Structural Plant Properties Using Radiative Transfer Models,
                                          <i>IGARSS 2024 - 2024 IEEE International Geoscience and Remote Sensing Symposium</i>, Athens, Greece, 2024, pp. 2811-2814,
                                          <a href="https://ieeexplore.ieee.org/document/10642442" target="_blank">doi: 10.1109/IGARSS53475.2024.10642442</a>.')),
                                         
                                         p(style = "text-align: justify;",HTML('Arano et al., (2024). Enhancing Chlorophyll Content Estimation With Sentinel-2 Imagery: A Fusion of Deep Learning and Biophysical Models,
                                          <i>IGARSS 2024 - 2024 IEEE International Geoscience and Remote Sensing Symposium</i>, Athens, Greece, 2024, pp. 4486-4489,
                                          <a href="https://ieeexplore.ieee.org/document/10641613" target="_blank">doi: 10.1109/IGARSS53475.2024.10641613</a>.')),
                                         
                                         p(style = "text-align: justify;",'Camino et al., (in prep). Integrating physiological plant traits with
                                             Sentinel-2 imagery for monitoring gross primary production and detecting forest disturbances. '),
                                         
                                         
                                ), #end TabPanel
                                
                                tabPanel("References for RTMs",
                                         
                                         br(),
                                         "If you utilize any of the RTM models within the ToolsRTM package, please be sure to cite the following references:",
                                         
                                         h4("PROSPECT model"),
                                         
                                         p(style = "text-align: justify;","Jacquemoud, S., Baret, F., 1990. PROSPECT: a model of leaf optical properties spectra. Remote Sens. Environ. 34, 75–91. ",
                                           br(),
                                           tags$a(href="https://doi.org/10.1016/0034-4257 (90)90100-Z", "https://doi.org/10.1016/0034-4257 (90)90100-Z")),
                                         p(style = "text-align: justify;","Jacquemoud S, Baret F, Hanocq J-F, 1992. Modeling spectral and bidirectional soil reflectance. Remote Sensing of Environment, 41, 123–132.",
                                           br(),
                                           tags$a(href="https://doi.org/10.1016/0034-4257(92)90072-R", "https://doi.org/10.1016/0034-4257(92)90072-R")),
                                         
                                         p(style = "text-align: justify;","Féret J-B, Gitelson AA, Noble SD & Jacquemoud S, 2017. PROSPECT-D: Towards modeling leaf optical properties through a complete lifecycle. Remote Sensing of Environment, 193, 204–215.",
                                           tags$a(href="https://doi.org/10.1016/j.rse.2017.03.004", "https://doi.org/10.1016/j.rse.2017.03.004")),
                                         p(style = "text-align: justify;","Féret, J.B., Berger, K., de Boissieu, F., Malenovský, Z., 2021. PROSPECT-PRO for estimating content of nitrogen-containing leaf proteins and other carbon-based constituents. Remote Sens. Environ. 252.",
                                           tags$a(href="https://doi.org/10.1016/j.rse.2020.112173", "https://doi.org/10.1016/j.rse.2020.112173")),
                                         
                                         h4("Liberty model"),
                                         
                                         p(style = "text-align: justify;","Dawson, T. P., Curran, P. J., & Plummer, S. E. (1998). LIBERTY—Modeling the Effects of Leaf Biochemical Concentration on Reflectance Spectra. Remote Sensing of Environment, 65(1), 50–60.",
                                           tags$a(href="https://doi.org/10.1016/S0034-4257(98)00007-8", "https://doi.org/10.1016/S0034-4257(98)00007-8")),
                                         p(style = "text-align: justify;","Di Vittorio, A. V. (2009). Enhancing a leaf radiative transfer model to estimate concentrations and in vivo specific absorption coefficients of total carotenoids and chlorophylls a and b from single-needle reflectance and transmittance. Remote Sensing of Environment, 113(9), 1948–1966.",
                                           tags$a(href="https://doi.org/10.1016/j.rse.2009.05.002", "https://doi.org/10.1016/j.rse.2009.05.002")),
                                         
                                         h4("FLUSPECT model"),
                                         
                                         p(style = "text-align: justify;","Vilfan, N., van der Tol, C., Muller, O., Rascher, U., Verhoef, W., 2016. Fluspect-B: A model for leaf fluorescence, reflectance and transmittancespectra.
                                           Remote Sens. Environ. 186, 596?615.",
                                           tags$a(href="https://doi:10.1016/j.rse.2016.09.017", "https://doi:10.1016/j.rse.2016.09.017")),
                                         
                                         h4("fourSAIL & fourSAIL-2 models"),
                                         
                                         p(style = "text-align: justify;","Verhoef W & Bach H, 2007. Coupled soil–leaf-canopy and atmosphere radiative transfer modeling to simulate hyperspectral multi-angular surface reflectance and TOA radiance data. Remote Sensing of Environment, 109:166-182.",
                                           tags$a(href="https://doi:10.1016/j.rse.2006.12.013", "https://doi:10.1016/j.rse.2006.12.013")),
                                         
                                         p(style = "text-align: justify;","Verhoef W, Jia L, Xiao Q & Su Z, 2007. Unified optical-thermal four-stream radiative transfer theory for homogeneous vegetation canopies. IEEE Transactions in Geosciences and Remote Sensing, 45:1808–1822.",
                                           tags$a(href=" https://doi.org/10.1109/TGRS.2007.895844", " https://doi.org/10.1109/TGRS.2007.895844")),
                                         
                                         p(style = "text-align: justify;","Jacquemoud S, Verhoef W, Baret F, Bacour C, Zarco-Tejada PJ, Asner GP, François C & Ustin SL, 2009. PROSPECT+ SAIL models: A review of use for vegetation characterization. Remote Sensing of Environment, 113:S56–S66. ",
                                           tags$a(href="https://doi.org/doi:10.1016/j.rse.2008.01.026", "https://doi.org/doi:10.1016/j.rse.2008.01.026")),
                                         
                                         h4("Invertible Forest Reflectance  Model"),
                                         
                                         p(style = "text-align: justify;",'Atzberger, C., 2000. Development of an Invertible Forest Reflectance Model: The INFOR- model.'),
                                         
                                         tags$p("Atzberger, C. (2000). Development of an invertible forest reflectance model: The INFOR-Model. In Buchroithner (Ed.), A decade of trans-european remote sensing cooperation. Proceedings of the 20th EARSeL Symposium Dresden, Germany, 14.-16. June 2000 (pp. 39-44)."),
                                         #br(),
                                         p(style = "text-align: justify;","Schlerf, M., Atzberger, C., 2006. Inversion of a forest reflectance model to estimate structural canopy variables from hyperspectral remote sensing data. Remote Sens. Environ. 100, 281–294. ",
                                           tags$a(href="https://doi.org/10.1016/j.rse.2005.10.006", "https://doi.org/10.1016/j.rse.2005.10.006")),
                                         #br(),
                                         #p("Schlerf, M., Atzberger, C., 2006. Inversion of a forest reflectance model to estimate structural canopy variables from hyperspectral remote sensing data. Remote Sens. Environ. 100, 281–294.",
                                         # tags$a(href='https://doi.org/10.1016/j.rse.2005.10.006.','https://doi.org/10.1016/j.rse.2005.10.006.')),
                                         #p(" ",
                                         #tags$a(href="", "")),
                                         br(),
                                         
                                ),# end Tab Reference panel
                                
                                tabPanel("About ToolsRTM", "",
                                         br(),
                                         p(style = "text-align: justify;",HTML('In this online RT-Simulator, we integrate the <strong>ToolsRTM</strong> package to simulate canopy reflectance using the primary radiative transfer (RT) models.
                                         This package allows for the rescaling of spectral resolution to accommodate various scales, including hyperspectral, Landsat and Sentinel-2.')),
                                         br(),
                                         
                                         p(style = "text-align: justify;",HTML('The <strong>ToolsRTM</strong> package is a vital component of the ToolsRTM Simulator app, which relies on two key R packages to ensure optimal functionality.
                                         These packages are essential for the efficient simulation and execution of functions within the server.')),
                                         br(),
                                         p(style = "text-align: justify;",HTML('<strong>ToolsRTM</strong> Package: This comprehensive package consolidates a variety of RT models for simulating reflectance at the Top of Canopy (TOC).
                                         It features several canopy models, including fourSAIL, fourSAIL2, and INFORM, as well as various leaf models. Additionally, the <strong>ToolsRTM</strong> package incorporates the Soil-Plant-Atmosphere Radiative Transfer (SPART) model, which simulates Top-of-Atmosphere (TOA) reflectance by integrating atmospheric parameters.
                                         SPART utilizes three computationally efficient sub-models—BSM for soil, PROSAIL for vegetation canopies, and SMAC for atmospheric effects—ensuring
                                         accurate simulation of directional TOA observations.')),
                                         br(),
                                         "Manual is available at ",
                                         tags$a(href="https://carlos-camino.shinyapps.io/0-toolsrtm-simulator/_w_e851c6b2/Notebooks/R/ToolsRTM/ToolsRTM.html",
                                                "ReadTheDocs"),
                                         ".",
                                         br(),
                                         br(),
                                         tags$figure(
                                           tags$img(src = "images/rtm_sims.png",height = "450px", width = "500px"),
                                           tags$figcaption("", style = "font-weight: bold;")
                                           
                                         ),
                                         strong("Fig 1."),'Reflectance simulation using several RT models.' ,
                                         
                                         #tags$img(src = "images/rtm_sims.png", height = "450px", width = "500px"),
                                         #htmltools::img(src = "images/rtm_sims.png"),
                                         
                                         br(),
                                         h4("Install ToolsRTM"),
                                         p(style = "text-align: justify;",HTML('<strong>ToolsRTM</strong> is avalaible on gitlab, so you can install using the R console:')),
                                         
                                         code('install.packages("toolsrtm-main.tar.gz",repos = NULL,type = "source")'),
                                         
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
                                         br(),
                                         
                                         
                                         br(),
                                         br(),
                                         br(),
                                         br(),
                                         
                                         #textInput("Install_package", "", value = "...."),
                                         
                                ), # end Tab Reference panel
                                
                              ) # close tabsetPanel
                              
                            ) # Close mainPanel
                          ) # Close sideBarLayout
                 ), # Close the Interactive ToolsRTM tab panel
)
# Define server logic required to draw a histogram ----
server <- function(input, output,session) {
  
  
  # 1) ToolsRTM   ---------------------------------------------------
  
  
  # Reactive expression for the leaf and canopy parameters -------------------------------
  
  params <- reactive({
    leaf_params <- switch(input$leaf_model,
                          
                          "PROSPECT-PRO" = c(N=input$N,Cab=input$Cab, Car=input$Car, Anth=input$Anth, Cbrown=input$Cbrown,
                                             EWT=input$EWT, LMA=0, alpha=40, #input$alpha,
                                             Prot=input$Prot, CBC=input$CBC),
                          "PROSPECT-D" = c(N=input$N_d,Cab=input$Cab_d, Car=input$Car_d,  Anth=input$Anth_d,Cbrown=input$Cbrown_d,
                                           EWT=input$EWT_d, LMA=input$LMA_d, alpha=40),#input$alpha_d),
                          
                          "Liberty" = c(cell.d=input$cell_d, inter.c=input$inter_c,
                                        baseline.abs=input$baseline_abs, leaf.thick=input$leaf_thick,
                                        albino.abs=input$albino_abs,
                                        Cab=input$Cab_l, EWT=input$EWT_l,
                                        lign.cell=input$lign_cell, Nitrogen=input$Nitrogen),
                          
                          
                          "FLUSPECT-B" = c(N=input$N_fd,Cab = input$Cab_fd, Car = input$Car_fd, Cs=input$Cs_fd,
                                           EWT=input$EWT_fd,LMA=input$LMA_fd,
                                           Cx = input$Cx_fd,
                                           fqe = input$fqe_fd),
                          
                          "FLUSPECT-B-Cx" = c(N=input$N_fp,Cab = input$Cab_fp, Car = input$Car_fp, Cs=input$Cs_fp,
                                              EWT=input$EWT_fp,LMA=input$LMA_fp,
                                              Cx = input$Cx_fp,
                                              fqe = input$fqe_fp))
    
    
    ## parameters for fourSAIL
    
    
    canopy_params <- switch(input$canopy_model,
                            "fourSAILH" = c(LAI=input$LAI, TypeLidf=2, LIDFa=input$LIDFa, LIDFb=0, hspot=input$hotspot,
                                            tts=input$tts, tto=input$tto, psi=input$psi, psoil=input$psoil),
                            
                            
                            "INFORM" = c(LAI=input$LAI_, TypeLidf=2, LIDFa=input$LIDFa_, LIDFb=0, hspot=input$hotspot_,
                                         tts=input$tts_, tto=input$tto_, psi=input$psi_, psoil=input$psoil_,
                                         LAIu=input$LAIu_, cd=input$cd_, sd=input$sd_, h=input$h_, skyl=0.1),
                            
                            
                            "fourSAILH2" = c(p1=input$p1_s2, p2=input$p2_s2, p3=input$p3_s2))
    return(list(leaf_params, canopy_params))
  })
  
  
  # Reactive expression for the LUT data ---------------------------------------------------
  lut_data.leaf <- reactive({
    
    # Get the selected leaf and canopy parameters from the reactive `params` expression
    leaf_params <- unlist(params()[[1]])
    leaf_params_names<-names(params()[[1]])
    
    # Call the LUT function using the selected parameters
    lut_leaf <- data.frame(rbind(leaf_params)) #, canopy_params))
    
    return(lut_leaf)
  })
  
  # Reactive expression for the LUT data ---------------------------------------------------
  lut_data.canopy <- reactive({
    # Get the selected leaf and canopy parameters from the reactive `params` expression
    canopy_params <- unlist(params()[[2]])
    # Call the LUT function using the selected parameters
    lut_canopy <- data.frame(rbind(canopy_params)) #, canopy_params))
    return(lut_canopy)
  })
  
  # Define the reactive expression for the leaf parameters ---------------------------------------------------
  leaf_params_prospect <- reactive({
    
    if(input$leaf_model == "PROSPECT-PRO" && (input$CBC != 0 & input$Prot != 0)) {
      
      leaf_params <- c(N = input$N, Cab = input$Cab, Car = input$Car, Anth = input$Anth,
                       Cbrown = input$Cbrown,EWT = input$EWT, LMA = 0, alpha = 40,#input$alpha,
                       Prot = input$Prot, CBC = input$CBC)
      
    } else if (input$leaf_model == "PROSPECT-D") {
      
      leaf_params <- c(N = input$N_d, Cab = input$Cab_d, Car = input$Car_d,  Anth=input$Anth_d, Cbrown = input$Cbrown_d,
                       alpha = 40, EWT = input$EWT_d, LMA = input$LMA_d) #input$alpha_d
      
    } else if (input$leaf_model == "Liberty") {
      
      leaf_params <- c(cell.d=input$cell_d, inter.c=input$inter_c,
                       baseline.abs=input$baseline_abs, leaf.thick=input$leaf_thick,
                       albino.abs=input$albino_abs,
                       Cab=input$Cab_l, EWT=input$EWT_l,
                       ## Adding PROSAIL-inputs for  Infinitive crown reflectance /unerstory
                       N = input$N, Car = input$Car, Anth = input$Anth,
                       Cbrown = input$Cbrown, LMA = input$LMA, alpha = 40,#input$alpha,
                       ##
                       lign.cell=input$lign_cell, Nitrogen=input$Nitrogen)
      
      
      
    } else if (input$leaf_model == "FLUSPECT-B") {
      
      leaf_params <- c(N = input$N_fd, Cab = input$Cab_fd, Car = input$Car_fd, Anth = 0, Cbrown= input$Cs_fd,
                       Cs = input$Cs_fd, Cx = input$Cx_fd,
                       EWT=input$EWT_fd,LMA=input$LMA_fd, alpha = 40, Prot = 0, CBC = 0,
                       fqe = input$fqe_fd)
      
      
      
    } else if (input$leaf_model == "FLUSPECT-B-Cx") {
      
      leaf_params <- c(N = input$N_fp, Cab = input$Cab_fp, Car = input$Car_fp, Anth=input$Anth_fp, Cbrown= input$Cs_fp,
                       Cs = input$Cs_fp, Cx = input$Cx_fp,
                       EWT=input$EWT_fp,LMA=input$LMA_fp, alpha = 40, Prot = input$Prot_fp,
                       CBC = input$CBC_fp,
                       fqe = input$fqe_fp)
      
      
      
    } else if (input$leaf_model == "PROSPECT-PRO" && (input$CBC == 0 & input$Prot == 0)) {
      
      leaf_params <- c(N = input$N, Cab = input$Cab, Car = input$Car, Anth = input$Anth, Cbrown = input$Cbrown,
                       EWT = input$EWT, LMA = input$LMA, alpha = 40, Prot = 0, CBC = 0) #input$alpha,
    }
    return(list(leaf_params))
    
  })
  
  
  # Reactive expression for the LUT data ---------------------------------------------------
  lut_data.sim <- reactive({
    # Get the selected leaf and canopy parameters from the reactive `params` expression
    leaf_params <- unlist(leaf_params_prospect()[[1]])
    
    canopy_params <- unlist(params()[[2]])
    
    # Call the LUT function using the selected parameters
    lut_leaf <- data.frame(rbind(leaf_params)) #, canopy_params))
    lut_canopy <- data.frame(rbind(canopy_params)) #, canopy_params))
    lut.to_sim <-cbind(lut_leaf,lut_canopy)
    
    
    # Call the Simulations
    data <- dataSpec_PDB
    Rsoil.dry  <- data[,11]  # rsoil1 = dry soil
    Rsoil.wet <- data[,12]  # rsoil2 = wet soil
    psoil = lut.to_sim[1,'psoil']
    #psoil	 <-  1    # soil factor (psoil=0: wet soil / psoil=1: dry soil)
    rsoil<- c(psoil*Rsoil.dry+(1-psoil)*Rsoil.wet)
    #reflectance_values <- -foursail(inputLUT=lut.to_sim[1,],rsoil=rsoil,LeafModel = 'PRO')
    return(list(lut.to_sim,rsoil))
    
  })
  
  
  
  # Reactive expression for the reflectance data ---------------------------------------------------
  reflectance_data <- reactive({
    
    LUT_ <- lut_data.sim()[[1]]
    print(LUT_)
    rsoil_ <- lut_data.sim()[[2]]
    if (input$canopy_model == 'fourSAILH'){
      LUT_$psi <- input$psi
      LUT_$psoil <- input$psoil
    } else {
      LUT_$psi <- input$psi_
      LUT_$psoil <- input$psoil
    }
    data <- dataSpec_PDB
    Rsoil.dry  <- data[,11]  # rsoil1 = dry soil
    Rsoil.wet <- data[,12]
    psoil = LUT_[1,'psoil']
    #psoil	 <-  1    # soil factor (psoil=0: wet soil / psoil=1: dry soil)
    rsoil_<- c(psoil*Rsoil.dry+(1-psoil)*Rsoil.wet)
    
    # Call the RTM function using the selected leaf and canopy parameters
    # return a list of two vectors
    if (input$canopy_model == 'fourSAILH' & input$leaf_model == 'PROSPECT-PRO'){
      
      reflectance_values <- foursail(inputLUT=LUT_[1,],rsoil=rsoil_, LeafModel = 'PROSPECT-PRO')
      rdot<-reflectance_values[[1]]
      rsot<-reflectance_values[[2]]
      # Compute the Bidirectional Reflectance Factor (BRF)
      reflectance_values<- Compute_BRF(rdot=rdot,rsot=rsot,tts=LUT_[1,'tts'],data.light=dataSpec_PDB)
      showNotification("Forward simulation with BRDF model successfully applied.", type = "message")
      
    } else if  (input$canopy_model == 'fourSAILH' & input$leaf_model == 'PROSPECT-D'){
      
      reflectance_values <- foursail(inputLUT=LUT_[1,],rsoil=rsoil_, LeafModel = 'PROSPECT-D')
      rdot<-reflectance_values[[1]]
      rsot<-reflectance_values[[2]]
      reflectance_values<- Compute_BRF(rdot=rdot,rsot=rsot,tts=LUT_[1,'tts'],data.light=dataSpec_PDB)
      showNotification("Forward simulation with BRDF model successfully applied.", type = "message")
      
    } else if  (input$canopy_model == 'fourSAILH' & input$leaf_model == 'Liberty'){
      
      reflectance_values <- foursail(inputLUT=LUT_[1,],rsoil=rsoil_,LeafModel = 'Liberty')
      rdot<-reflectance_values[[1]]
      rsot<-reflectance_values[[2]]
      reflectance_values<- Compute_BRF(rdot=rdot,rsot=rsot,tts=LUT_[1,'tts'],data.light=dataSpec_PDB)
      showNotification("Forward simulation with BRDF model successfully applied.", type = "message")
      
    } else if  (input$canopy_model == 'fourSAILH' & input$leaf_model == 'FLUSPECT-B'){
      
      reflectance_values <- foursail(inputLUT=LUT_[1,],rsoil=rsoil_[1:2001],LeafModel = 'Fluspect-B')
      rdot<-reflectance_values[[1]]
      rsot<-reflectance_values[[2]]
      reflectance_values<- Compute_BRF(rdot=rdot,rsot=rsot,tts=LUT_[1,'tts'],data.light=dataSpec_PDB, short.waves = T)
      showNotification("Forward simulation with BRDF model successfully applied.", type = "message")
      
    } else if  (input$canopy_model == 'fourSAILH' & input$leaf_model == 'FLUSPECT-B-Cx'){
      
      reflectance_values <- foursail(inputLUT=LUT_[1,],rsoil=rsoil_[1:2001],LeafModel = 'Fluspect-B-Cx')
      rdot<-reflectance_values[[1]]
      rsot<-reflectance_values[[2]]
      reflectance_values<- Compute_BRF(rdot=rdot,rsot=rsot,tts=LUT_[1,'tts'], short.waves = T)
      showNotification("Forward simulation with BRDF model successfully applied.", type = "message")
      
    } else if  (input$canopy_model == 'INFORM' & input$leaf_model == 'PROSPECT-PRO'){
      
      reflectance_values <- inform(inputLUT = LUT_[1,],rsoil=rsoil_,LeafModel = 'PROSPECT-PRO')
      showNotification("Forward simulation with INFORM model successfully applied.", type = "message")
      
      
    } else if  (input$canopy_model == 'INFORM' & input$leaf_model == 'PROSPECT-D'){
      
      reflectance_values <- inform(inputLUT = LUT_[1,],rsoil=rsoil_,LeafModel = 'PROSPECT-D')
      showNotification("Forward simulation with INFORM model successfully applied.", type = "message")
      
    } else if  (input$canopy_model == 'INFORM' & input$leaf_model == 'FLUSPECT-B'){
      
      reflectance_values <- inform(inputLUT = LUT_[1,],rsoil=rsoil_[1:2001],LeafModel = 'Fluspect-B')
      showNotification("Forward simulation with INFORM model successfully applied.", type = "message")
      
    } else if  (input$canopy_model == 'INFORM' & input$leaf_model == 'FLUSPECT-B-Cx'){
      
      reflectance_values <- inform(inputLUT = LUT_[1,],rsoil=rsoil_[1:2001],LeafModel = 'Fluspect-B-Cx')
      showNotification("Forward simulation with INFORM model successfully applied.", type = "message")
      
    } else if  (input$canopy_model == 'INFORM' & input$leaf_model == 'Liberty'){
      
      reflectance_values <- inform(inputLUT = LUT_[1,],rsoil=rsoil_,LeafModel = 'Liberty')
      showNotification("Forward simulation with INFORM model successfully applied.", type = "message")
      
    } else {
      
      reflectance_values<-sim.liberty$refl
      
    }
    wavelength=seq(400,2500,1)
    
    # Combine the reflectance values with wavelength values
    if (input$leaf_model == 'FLUSPECT-B-Cx' | input$leaf_model == 'FLUSPECT-B'){
      
      reflectance_df <- data.frame(wavelength = wavelength[1:2001], reflectance = reflectance_values)
      return(reflectance_df)
    } else {
      
      reflectance_df <- data.frame(wavelength = wavelength, reflectance = reflectance_values)
      return(reflectance_df)
    }
    
    
    
  })
  
  
  
  
  
  # Render the reflectance plot ---------------------------------------------------
  output$reflectance_plot <- renderPlot({
    
    
    
    if (input$sensor == 'RTM'){
      
      to_plot<-reflectance_data()
      
      
      #  if (input$checkbox_spectra) {
      
      #   to_plot_nd<-c(reflectance_data()$reflectance)
      
      
      #}
      
    } else if (input$sensor == 'Sentinel2a') {
      
      to_plot <- reflectance_data()
      
      
      Spec.simula<- data.frame(wave = to_plot$wavelength, rfl=to_plot$reflectance)
      sensor.i = Sentinel2A.MSI
      df_resampled <-get.spectral.convolution.rfl(df = Spec.simula,sensor.i, get.plots=F)
      to_plot <- data.frame(wavelength = df_resampled$wave, reflectance = df_resampled$RFL)
      
      showNotification("Satellite selection with applied bandwidths completed successfully.", type = "message")
      
      
      
      
      
    }  else if (input$sensor == 'Sentinel2b') {
      
      to_plot <- reflectance_data()
      
      Spec.simula<- data.frame(wave = to_plot$wavelength, rfl=to_plot$reflectance)
      sensor.i = Sentinel2B.MSI
      df_resampled <-get.spectral.convolution.rfl(df = Spec.simula,sensor.i, get.plots=F)
      to_plot <- data.frame(wavelength = df_resampled$wave, reflectance = df_resampled$RFL)
      showNotification("Satellite selection with applied bandwidths completed successfully.", type = "message")
      
      
    }   else if (input$sensor == 'Landsat4') {
      
      to_plot <- reflectance_data()
      
      
      Spec.simula<- data.frame(wave = to_plot$wavelength, rfl=to_plot$reflectance)
      sensor.i = LANDSAT4.TM
      df_resampled <-get.spectral.convolution.rfl(df = Spec.simula,sensor.i, get.plots=F)
      to_plot <- data.frame(wavelength = df_resampled$wave, reflectance = df_resampled$RFL)
      showNotification("Satellite selection with applied bandwidths completed successfully.", type = "message")
      
      
    } else if (input$sensor == 'Landsat5') {
      
      to_plot <- reflectance_data()
      
      Spec.simula<- data.frame(wave = to_plot$wavelength, rfl=to_plot$reflectance)
      sensor.i = LANDSAT5.TM
      df_resampled <-get.spectral.convolution.rfl(df = Spec.simula,sensor.i, get.plots=F)
      to_plot <- data.frame(wavelength = df_resampled$wave, reflectance = df_resampled$RFL)
      showNotification("Satellite selection with applied bandwidths completed successfully.", type = "message")
      
      
      
    } else if (input$sensor == 'Landsat7') {
      
      to_plot <- reflectance_data()
      
      Spec.simula<- data.frame(wave = to_plot$wavelength, rfl=to_plot$reflectance)
      sensor.i = LANDSAT7.ETM
      df_resampled <-get.spectral.convolution.rfl(df = Spec.simula,sensor.i, get.plots=F)
      to_plot <- data.frame(wavelength = df_resampled$wave, reflectance = df_resampled$RFL)
      showNotification("Satellite selection with applied bandwidths completed successfully.", type = "message")
      
      
    } else if (input$sensor == 'Landsat8') {
      
      to_plot <- reflectance_data()
      
      Spec.simula<- data.frame(wave = to_plot$wavelength, rfl=to_plot$reflectance)
      sensor.i = LANDSAT8.OLI
      df_resampled <-get.spectral.convolution.rfl(df = Spec.simula,sensor.i, get.plots=F)
      to_plot <- data.frame(wavelength = df_resampled$wave, reflectance = df_resampled$RFL)
      showNotification("Satellite selection with applied bandwidths completed successfully.", type = "message")
      
      
    } else if (input$sensor == 'MODIS') {
      
      to_plot <- reflectance_data()
      showNotification("Satellite selection displayed native resolution at 1nm.", type = "warning")
      showNotification("Spectral resampling for this sensor is currently under development.", type = "error")
      
    } else if (input$sensor == 'RapidEye') {
      
      to_plot <- reflectance_data()
      showNotification("Satellite selection displayed native resolution at 1nm.", type = "warning")
      showNotification("Spectral resampling for this sensor is currently under development.", type = "error")
      
    } else if (input$sensor == 'Quickbird') {
      
      to_plot <- reflectance_data()
      showNotification("Satellite selection displayed native resolution at 1nm.", type = "warning")
      showNotification("Spectral resampling for this sensor is currently under development.", type = "error")
      
    } else if (input$sensor == 'EnMAP') {
      
      to_plot <- reflectance_data()
      showNotification("Satellite selection displayed native resolution at 1nm.", type = "warning")
      showNotification("Spectral resampling for this sensor is currently under development.", type = "error")
      
    } else if (input$sensor == 'PRISMA') {
      
      to_plot <- reflectance_data()
      showNotification("Satellite selection displayed native resolution at 1nm.", type = "warning")
      showNotification("Spectral resampling for this sensor is currently under development.", type = "error")
      
    } else if (input$sensor == 'Hyperion') {
      
      to_plot <- reflectance_data()
      showNotification("Satellite selection displayed native resolution at 1nm.", type = "warning")
      showNotification("Spectral resampling for this sensor is currently under development.", type = "error")
      
    } else if (input$sensor == 'WorldView2-4') {
      
      to_plot <- reflectance_data()
      showNotification("Satellite selection displayed native resolution at 1nm.", type = "warning")
      showNotification("Spectral resampling for this sensor is currently under development.", type = "error")
      
      
    } else if (input$sensor == 'WorldView2-8') {
      
      to_plot <- reflectance_data()
      showNotification("Satellite selection displayed native resolution at 1nm.", type = "warning")
      showNotification("Spectral resampling for this sensor is currently under development.", type = "error")
      
    } else if (input$sensor == 'ALI') {
      
      to_plot <- reflectance_data()
      showNotification("Satellite selection displayed native resolution at 1nm.", type = "warning")
      showNotification("Spectral resampling for this sensor is currently under development.", type = "error")
      
      
      
    }
    
    # Plot Polygons by Sensor  ---------------------------------------------------
    
    # Create a list of polygons for each bandset ---------------------------------------------------
    sentinel2a_polygons <- list(
      geom_rect(xmin = 432.2, xmax = 453.2, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 459.4, xmax = 525.4, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 541.8, xmax = 577.8, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 649.1, xmax = 680.1, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 696.6, xmax = 711.6, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 733.0, xmax = 748.0, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 772.8, xmax = 792.8, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 779.8, xmax = 885.8, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 854.2, xmax = 875.2, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 935.1, xmax =  955.1, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 1358.0, xmax = 1389.0, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 1568.2, xmax =   1659.2, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 2114.9, xmax =  2289.9, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5))
    
    sentinel2b_polygons <- list(
      geom_rect(xmin = 431.8, xmax = 452.8, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 459.1, xmax = 525.1, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 541.0, xmax = 577.0, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 649.5, xmax = 680.5, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 695.8, xmax = 711.8, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 731.6, xmax = 746.6, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 769.7, xmax = 789.7, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 780.0, xmax = 886.0, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 853.0, xmax = 875.0, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 932.7, xmax =  953.7, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 1361.9, xmax = 1391.9, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 1563.4, xmax = 1657.4, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 2093.2, xmax =  2278.2, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5))
    
    modis_polygons <- list(
      geom_rect(xmin = 405, xmax = 420, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 438, xmax = 448, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 459, xmax = 479, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 483, xmax = 493, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 526, xmax = 536, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 546, xmax = 556, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 545, xmax = 565, ymin = -Inf, ymax = Inf, fill = "grey90", alpha = 0.5),
      geom_rect(xmin = 662, xmax = 672, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 673, xmax = 683, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 743, xmax = 753, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 862, xmax = 877, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 890, xmax =  920, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 915, xmax = 965, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 931, xmax = 941, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 1230, xmax = 1250, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 1628, xmax = 1652, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 2105, xmax = 2155, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5)
    )
    
    landsat45_polygons <- list(
      geom_rect(xmin = 450, xmax = 520, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 520, xmax = 600, ymin = -Inf, ymax = Inf, fill = "grey90", alpha = 0.5),
      geom_rect(xmin = 630, xmax = 690, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 760, xmax = 900, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 1550, xmax = 1750, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 2080, xmax = 2350, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5))
    
    landsat7_polygons <- list(
      geom_rect(xmin = 450, xmax = 520, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 530, xmax = 610, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 630, xmax = 690, ymin = -Inf, ymax = Inf, fill = "grey90", alpha = 0.5),
      geom_rect(xmin = 780, xmax = 900, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 1550, xmax = 1750, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 2090, xmax = 2350, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5))
    
    landsat8_polygons <- list(
      geom_rect(xmin = 427, xmax = 462, ymin = -Inf, ymax = Inf, fill = "grey90", alpha = 0.5),
      geom_rect(xmin = 435, xmax = 530, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 512, xmax = 602, ymin = -Inf, ymax = Inf, fill = "grey90", alpha = 0.5),
      geom_rect(xmin = 625, xmax = 685, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 830, xmax = 897, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 1340, xmax = 1405, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 1515, xmax = 1697, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 2037, xmax = 2352, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5))
    
    ali_polygons <- list(
      geom_rect(xmin = 433, xmax = 453, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 450, xmax = 525, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 525, xmax = 605, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 630, xmax = 690, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 775, xmax = 805, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 845, xmax = 890, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 1200, xmax = 1300, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 1550, xmax = 1750, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5),
      geom_rect(xmin = 2080, xmax = 2350, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5))
    
    if (input$sensor != 'RTM') {
      
      
      
      # Is Add region is active, ggplot plot the spectral bandwidth  of the sensors
      
      if (input$AddRegion == T) {
        
        showNotification("Please note that this process may require additional time.", type = "warning")
        
        # Determine which polygons to add based on the selected sensor
        if (input$sensor == "Sentinel2a") {
          polygons_to_add <- sentinel2a_polygons
        } else if (input$sensor == "Sentinel2b") {
          polygons_to_add <- sentinel2b_polygons
        } else if (input$sensor == "Landsat4" | input$sensor == "Landsat5") {
          polygons_to_add <- landsat45_polygons
        } else if (input$sensor == "Landsat7") {
          polygons_to_add <- landsat7_polygons
        } else if (input$sensor == "Landsat8") {
          polygons_to_add <- landsat8_polygons
        } else if (input$sensor == "ALI") {
          polygons_to_add <- ali_polygons
        } else if (input$sensor == "MODIS") {
          polygons_to_add <- modis_polygons
        } else {
          polygons_to_add <- NULL
        }
        
      } else {
        polygons_to_add <- NULL
      }
      
      
      plot_satellite <-ggplot(reflectance_data(), aes(x = wavelength, y = reflectance)) +
        polygons_to_add +
        geom_line(aes(color='RTM'), size = 0.8) +
        
        geom_line(data = to_plot, aes(x = wavelength, y = reflectance, color = "Sensor"),
                  lwd=1, linetype = "dashed") +
        geom_point(data = to_plot, aes(x = wavelength, y = reflectance), color = "red",cex=3) +
        scale_color_manual(name = "Reflectance", values = c("RTM" = "black", "Sensor" = "dodgerblue3"),
                           labels = c("RTM", input$sensor)) +
        labs(x = "Wavelength (nm)", y = "Reflectance") +
        # Set axis labels and title
        theme_bw() + ylim(0,0.7) +
        theme(legend.position = "insider",
              legend.box.background = element_rect(color = "black",linewidth=1),
              plot.title = element_text(hjust = 0.5, size=14,face="bold"),
              axis.title = element_text(face="bold", size=14),
              legend.text = element_text(face="bold", size=10),
              axis.text.y=element_text(hjust = 0.5, size=12,face="bold"),
              axis.text.x=element_text(hjust = 0.5, size=12,face="bold"),
              legend.title=element_blank()) +
        # legend.title = element_text(face = "bold", size = 14)) +
        guides(color = guide_legend(title.position = "top", title.hjust = 0.5))
      
      plot_satellite
      
      
    } else {
      
      plot_satellite <-ggplot(to_plot, aes(x = wavelength, y = reflectance)) +
        geom_line(aes(color='RTM'), linewidth = 1,color='black') +
        labs(x = "Wavelength (nm)", y = "Reflectance") +
        theme_bw() + ylim(0,0.7) +
        theme(legend.position = 'none',
              legend.box.background = element_rect(color = "black",linewidth=1),
              plot.title = element_text(hjust = 0.5, size=14,face="bold"),
              axis.title = element_text(face="bold", size=14),
              legend.text = element_text(face="bold", size=10),
              axis.text.y=element_text(hjust = 0.5, size=12,face="bold"),
              axis.text.x=element_text(hjust = 0.5, size=12,face="bold"),
              legend.title=element_blank()) +
        # legend.title = element_text(face = "bold", size = 14)) +
        guides(color = guide_legend(title.position = "top", title.hjust = 0.5))
      
      plot_satellite
      
    }
    
    
  })
  # Plot Polygons by Sensor  -------------------------
  
  
  # Create a reactive data frame for saving data  ---------------------------------------------------
  
  export_table <- reactive({
    
    if (input$sensor != 'RTM') {
      rfl_ <- t(reflectance_data()$reflectance)
      wave_names <- paste('R',reflectance_data()$wavelength,sep='.')
    } else {
      rfl_ <- t(reflectance_data()$reflectance)
      wave_names <- paste('R',reflectance_data()$wavelength,sep='.')
      
    }
    
    LUT_ <- lut_data.sim()[[1]]
    to_save<-cbind(LUT_,rfl_)
    # Assign column names if necessary
    colnames(to_save) <- c(names(LUT_), wave_names)  # Adjust as per your LUT column names
    return(to_save)
    
  })
  #### Add the download handler for saving data ---------------
  output$downloadData_sim <- downloadHandler(
    filename = function() {
      paste("LUT_", input$leaf_model,'_', input$canopy_model,'_', input$sensor, ".csv", sep = "")
    },
    content = function(file) {
      
      write.csv(export_table(), file, row.names = FALSE)
      showNotification("Simulation at 1nm saved successfully.", type = "message")
      
    }
  )
  
  # Inside your server function
  observeEvent(input$checkbox, {
    if (input$checkbox) {
      # showNotification("Reflectance spectra will be accumulated.", type = "message")
      showNotification("Reflectance spectra accumulation is currently disabled. We are working on enabling this feature.", type = "warning")
    }
  })
  
  
  # Inside your server function
  observeEvent(input$checkbox_sif, {
    if (input$checkbox_sif) {
      # showNotification("Reflectance spectra will be accumulated.", type = "message")
      showNotification("SIF spectra accumulation is currently disabled. We are working on enabling this feature.", type = "warning")
    }
  })
  
  # Inside your server function
  observeEvent(input$checkbox_scope, {
    if (input$checkbox_scope) {
      # showNotification("Reflectance spectra will be accumulated.", type = "message")
      showNotification("Reflectance spectra accumulation is currently disabled. We are working on enabling this feature.", type = "warning")
    }
  })
  # Inside your server function
  observeEvent(input$checkbox_LO, {
    if (input$checkbox_LO) {
      # showNotification("Reflectance spectra will be accumulated.", type = "message")
      showNotification("Radiance spectra accumulation is currently disabled. We are working on enabling this feature.", type = "warning")
    }
  })
  
  ### Define table for Liberty model ---------------------------------------------------
  
  output$Liberty_table <- renderTable({
    data <- data.frame(
      Parameter = c('cell.d','inter.c','baseline.abs', 'leaf.thick',
                    'albino_abs','Cab','EWT','lign.cell','Nitrogen'),
      Description = c('Cell diameter','Intercellular air space', 'baseline absorption','leaf thickness',
                      'Albino absorption', 'Chlorophyll content', 'Leaf water content','Lignin and cellulose content',
                      'leaf nitrogen'),
      Units = c('m-6','-','-','-','-','μg cm-2','μg cm-2','-','g cm-2'),
      Min = c(20,0.01,0.0004,1,0,0,0,10,0.3),
      Max = c(200,0.1,0.0006,10,4,60,0.05,80,2),
      Default = c(45,0.0045,0.0004,1.6,2,40,0.009,40,1)
    )
    colnames(data) <- c("Input", "Description", "Units", "Min", "Max","default")
    return(data)
  })
  
  
  ### Define table for Prospect model ---------------------------------------------------
  
  
  output$prospect_table <- renderTable({
    data <- data.frame(
      Parameter = c('N','Cab','Car', 'Anth','Cbrown','EWT','LMA','alpha','Prot','CBC'),
      Description = c('Leaf mesophyll','Chlorophyll content', 'Carotenoid content', 'Anthocyanins content',
                      'Brown pigments','Leaf water content','Leaf matter content','alpha','leaf proteins','carbon-based Constutient'),
      Units = c('-','μg cm-2','μg cm-2','μg cm-2','-','g cm-2','g cm-2','-','g cm-2','g cm-2'),
      Min = c(1,0,0,0,0,0,0,0,0,0),
      Max = c(4,100,40,7,1,0.05,0.05,60,0.05,0.05),
      Default = c(2.5,40,10,0,0,0.009,0.0012,40,0.01,0.01)
    )
    colnames(data) <- c("Input", "Description", "Units", "Min", "Max","default")
    return(data)
  })
  
  
  ### Define table for Prospect model ---------------------------------------------------
  
  
  output$fluspect_table <- renderTable({
    data <- data.frame(
      Parameter = c('N','Cab','Car', 'Anth','Cbrown','EWT','LMA','Cs','Prot','CBC', 'Cx'),
      Description = c('Leaf mesophyll','Chlorophyll content', 'Carotenoid content', 'Anthocyanins content',
                      'Brown pigments','Leaf water content','Leaf matter content','Leaf senescence',
                      'leaf proteins','carbon-based Constutient', 'Violaxanthin - Zeaxanthin transition status'),
      Units = c('-','μg cm-2','μg cm-2','μg cm-2','-','g cm-2','g cm-2','-','g cm-2','g cm-2','-'),
      Min = c(1,0,0,0,0,0,0,0,0,0,0),
      Max = c(4,100,40,7,1,0.05,0.05,1,0.05,0.05,1),
      Default = c(2.5,40,10,0,0,0.009,0.0012,0.1,0.01,0.01,0.1)
    )
    colnames(data) <- c("Input", "Description", "Units", "Min", "Max","default")
    return(data)
  })
  
  
  
  ### Define table for INFORM  model ---------------------------------------
  
  output$inform_table <- renderTable({
    data <- data.frame(
      Parameter = c('LAIu','sd','h', 'cd'),
      Description = c('Leaf Area Index understorey','Stem density ','Tree height',
                      'Crown diameter'),
      Units = c('m2 m-2','ha-1','m','m'),
      Min = c(0.01,0,0,0.1),
      Max = c(3,3000,50,10),
      Default = c(0.1,650,20,4.5)
    )
    colnames(data) <- c("Input", "Description", "Units", "Min", "Max","default")
    return(data)
  })
  
  
  
  ### Define table for foursailh model ---------------------------------------------------
  
  output$foursailh_table <- renderTable({
    data <- data.frame(
      Parameter = c('LAI','LIDFa','LIDFb', 'hotspot',
                    'tts','tto','psi','psoil','Soil reflectance'),
      Description = c('Leaf Area Index','leaf inclination distribution function a','Type of leaf inclination distribution function b',
                      'Hot Spot parameter','Sun zeith angle','Observer zeith angle','zimuth Sun / Observer','soil factor','Soil reflectance'),
      Units = c('m2 m-2','deg','deg','-','deg','deg','deg','-','%'),
      Min = c(0.001,0,0,0,0,0,0,0,NA),
      Max = c(10,90,0,1,90,90,180,1,NA),
      Default = c(4,30,0,0.5,0,30,0,0.2,NA)
    )
    colnames(data) <- c("Input", "Description", "Units", "Min", "Max","default")
    return(data)
  })
}

# Create Shiny object
shinyApp(ui = ui, server = server)

