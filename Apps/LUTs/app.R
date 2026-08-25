
rm(list= ls())

# 1. load the main libraries  -----

library(ToolsRTM)
required_packages <- c("shiny", "shinythemes", "shinybusy","ggplot2", "dplyr", 'tidyverse',
                       "doParallel",'parallel','foreach','DT')

# Check for missing packages and install them if necessary
missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}
# Load the libraries
lapply(required_packages, library, character.only = TRUE)

# Define UI

ui <- navbarPage("Online reflectance simulator",theme = shinytheme("flatly"),
                 
                 tabPanel("Create your Look-up table", # Open a second tabPanel for generating LUTs
                        
                          sidebarLayout(
                            sidebarPanel(
                              selectInput("leaf_lut", "Select a leaf RT Model:",
                                          choices = c("PROSPECT-PRO", "PROSPECT-D","Liberty","FLUSPECT-B","FLUSPECT-B-Cx")), # ,""
                              selectInput("canopy_lut", "Select a canopy RT Model:",
                                          choices = c('fourSAILH','INFORM')), ##fourSAILH2
                              selectInput("sensor_lut", "Select sensor for resampling resolution:",
                                          choices = c('RTM','Sentinel2a','Sentinel2b','Landsat4','Landsat5','Landsat7','Landsat8')),
                              numericInput("n_samples_lut", "Number of samples:", value = 25,min = 10,max = 20000),
                              p('maximum LUT to 20,000 simulations'),
                              numericInput("seed_lut", "ramdom seed:", value = 1234),
                              p('Random seed parameter for repeatability'),
                              h4(strrep("-", 60)),
                              h4('Leaf parameters:'),
                              
                              conditionalPanel(
                                condition = "input.leaf_lut == 'PROSPECT-PRO'",
                                sliderInput("Cab2", HTML("Chlorophyll content (μg cm<sup>-2</sup>)"), min = 0, max = 100, value = c(2,60)),
                                checkboxInput("Dist_cab2", "Use a gauss distribution", value = F),
                                
                                sliderInput("Car2", HTML("Carotenoid content (μg cm<sup>-2</sup>)"), min = 0, max = 40, value = c(5,20)),
                                checkboxInput("Dist_car2", "Use a gauss distribution", value = TRUE),
                                
                                sliderInput("Anth2", HTML("Anthocyanin content (μg cm<sup>-2</sup>)"), min = 0, max = 10, value = c(0,7)),
                                checkboxInput("Dist_anth2", "Use a gauss distribution", value = TRUE),
                                
                                sliderInput("Cbrown2", "Cbrown", min = 0, max = 1, value = c(0,1)),
                                checkboxInput("Dist_cbrown2", "Use a gauss distribution", value = FALSE),
                                
                                sliderInput("N2", "mesophyll structure parameter", min = 1, max = 4, value = c(1.5,3)),
                                checkboxInput("Dist_n2", "Use a gauss distribution", value = FALSE),
                                
                                sliderInput("EWT2", HTML("Water content (g cm<sup>-2</sup>)"), min = 0.0001, max = 0.05, value = c(0.01,0.03)),
                                checkboxInput("Dist_ewt2", "Use a gauss distribution", value = FALSE),
                                
                                sliderInput("Prot2", HTML("Proteins (g cm<sup>-2</sup>)"), min = 0.0001, max = 0.03, value = c(0.0001,0.005)),
                                checkboxInput("Dist_prot2", "Use a gauss distribution", value = FALSE),
                                
                                sliderInput("CBC2",  HTML("Carbon-based constituent (g cm<sup>-2</sup>)"), min = 0, max = 0.03, value = c(0,0.003)),
                                checkboxInput("Dist_cbc2", "Use a gauss distribution", value = FALSE),
                                
                              ),
                              
                              # Leaf Model parameters (4)
                              conditionalPanel(
                                condition = "input.leaf_lut == 'FLUSPECT-B'",
                                
                                sliderInput("fqe_fd_2", "Fluorescence quantum efficiency (fqe)", min = 0, max = 0.05, value = c(0.02,0.03)),
                                checkboxInput("Dist_fqe_fd_2", "Use a gauss distribution", value = T),
                                
                                sliderInput("Cab_fd_2", HTML("Chlorophyll content (μg cm<sup>-2</sup>)"), min = 0, max = 100,  value = c(10,50)),
                                checkboxInput("Dist_cab_fd_2", "Use a gauss distribution", value = F),
                                
                                sliderInput("Car_fd_2", HTML("Carotenoid content (μg cm<sup>-2</sup>)"), min = 0, max = 20,  value = c(4,10)),
                                checkboxInput("Dist_car_fd_2", "Use a gauss distribution", value = T),
                                
                                sliderInput("Cs_fd_2", "Leaf Senescence", min = 0, max = 1,  value = c(0.2,0.3)),
                                checkboxInput("Dist_cs_fd_2", "Use a gauss distribution", value = F),
                                
                                sliderInput("Cx_fd_2", "Violaxanthin - Zeaxanthin transition status", min = 0, max = 1,  value = c(0.5,0.8)),
                                checkboxInput("Dist_cx_fd_2", "Use a gauss distribution", value = F),
                                
                                sliderInput("N_fd_2", "mesophyll structure parameter", min = 1, max = 4,  value = c(1.5,2.5)),
                                checkboxInput("Dist_n_fd_2", "Use a gauss distribution", value = F),
                                
                                sliderInput("EWT_fd_2", HTML("Water content (g cm<sup>-2</sup>)"), min = 0, max = 0.05, value = c(0.01,0.02)),
                                checkboxInput("Dist_ewt_fd_2", "Use a gauss distribution", value = T),
                                
                                sliderInput("LMA_fd_2", HTML("dry matter content (g cm<sup>-2</sup>)"), min = 0, max = 0.05,  value = c(0.02,0.04)),
                                checkboxInput("Dist_lma_fd_2", "Use a gauss distribution", value = T),
                                
                                
                              ),
                              
                              # Leaf Model parameters (4)
                              conditionalPanel(
                                condition = "input.leaf_lut == 'FLUSPECT-B-Cx'",
                                
                                sliderInput("fqe_fp_2", "Fluorescence quantum efficiency (fqe)", min = 0, max = 0.05, value = c(0.01,0.02)),
                                checkboxInput("Dist_fqe_fp_2", "Use a gauss distribution", value = T),
                                
                                sliderInput("Cab_fp_2", HTML("Chlorophyll content (μg cm<sup>-2</sup>)"), min = 0, max = 100,  value = c(30,60)),
                                checkboxInput("Dist_cab_fp_2", "Use a gauss distribution", value = T),
                                
                                
                                sliderInput("Car_fp_2", HTML("Carotenoid content (μg cm<sup>-2</sup>)"), min = 0, max = 40,  value = c(10,20)),
                                checkboxInput("Dist_car_fp_2", "Use a gauss distribution", value = T),
                                
                                sliderInput("Anth_fp_2",  HTML("Anthocyanin content (μg cm<sup>-2</sup>)"), min = 0, max = 7,  value = c(0.5,6)),
                                checkboxInput("Dist_ant_fp_2", "Use a gauss distribution", value = T),
                                
                                sliderInput("Cs_fp_2", "Leaf Senescence", min = 0, max = 1,  value = c(0.1,0.5)),
                                checkboxInput("Dist_cs_fp_2", "Use a gauss distribution", value = T),
                                
                                sliderInput("Cx_fp_2", "Violaxanthin - Zeaxanthin transition status", min = 0, max = 1,  value = c(0.1,0.5)),
                                checkboxInput("Dist_cx_fp_2", "Use a gauss distribution", value = T),
                                
                                sliderInput("N_fp_2", "mesophyll structure parameter", min = 1, max = 4,  value = c(1.3,2.5)),
                                checkboxInput("Dist_n_fp_2", "Use a gauss distribution", value = T),
                                
                                sliderInput("EWT_fp_2", HTML("Water content (g cm<sup>-2</sup>)"), min = 0, max = 0.05,  value = c(0.01,0.02)),
                                checkboxInput("Dist_ewt_fp_2", "Use a gauss distribution", value = T),
                                
                                sliderInput("LMA_fp_2", HTML("dry matter content (g cm<sup>-2</sup>)"), min = 0, max = 0.05, value = c(0.01,0.02)),
                                checkboxInput("Dist_lma_fp_2", "Use a gauss distribution", value = T),
                                
                                sliderInput("Prot_fp_2", HTML("Proteins (g cm<sup>-2</sup>)"), min = 0, max = 0.03,  value = c(0.001,0.005)),
                                checkboxInput("Dist_prot_fp_2", "Use a gauss distribution", value = T),
                                
                                sliderInput("CBC_fp_2", HTML("Carbon-based Constituent  (g cm<sup>-2</sup>)"), min = 0, max = 0.03,  value = c(0.001,0.002)),
                                checkboxInput("Dist_cbc_fp_2", "Use a gauss distribution", value = T),
                                
                              ),
                              # Leaf Model parameters (2)
                              conditionalPanel(
                                condition = "input.leaf_lut == 'PROSPECT-D'",
                                
                                sliderInput("Cab_d2",  HTML("Chlorophyll content (μg cm<sup>-2</sup>)"), min = 0, max = 100, value = c(5,60)),
                                checkboxInput("Dist_cab_d2", "Use a gauss distribution", value = F),
                                
                                sliderInput("Car_d2", HTML("Carotenoid content (μg cm<sup>-2</sup>)"), min = 0, max = 40,  value = c(4,20)),
                                checkboxInput("Dist_car_d2", "Use a gauss distribution", value = F),
                                
                                sliderInput("Anth_d2",  HTML("Anthocyanin content (μg cm<sup>-2</sup>)"), min = 0, max = 10,  value = c(0,5)),
                                checkboxInput("Dist_anth_d2", "Use a gauss distribution", value = F),
                                
                                sliderInput("Cbrown_d2", "Cbrown", min = 0, max = 1,  value = c(0,1)),
                                checkboxInput("Dist_cbrown_d2", "Use a gauss distribution", value = F),
                                
                                sliderInput("N_d2", "mesophyll structure parameter", min = 1, max = 5,  value = c(1,3)),
                                checkboxInput("Dist_n_d2", "Use a gauss distribution", value = F),
                                
                                sliderInput("EWT_d2", HTML("Water content (g cm<sup>-2</sup>)"), min = 0, max = 0.05,  value = c(0.01,0.03)),
                                checkboxInput("Dist_ewt_d2", "Use a gauss distribution", value = F),
                                
                                sliderInput("LMA_d2", HTML("dry matter content (g cm<sup>-2</sup>)"), min = 0, max = 0.05,  value = c(0.01,0.03)),
                                checkboxInput("Dist_lma_d2", "Use a gauss distribution", value = F),
                                
                                
                                
                              ),
                              # Leaf Model parameters (3)
                              conditionalPanel(
                                condition = "input.leaf_lut == 'Liberty'",
                                
                                sliderInput("Cab_l2", HTML("Chlorophyll content (μg cm<sup>-2</sup>)"), min = 0, max = 80, value = c(40,70)),
                                checkboxInput("Dist_cab_l2", "Use a gauss distribution", value = F),
                                
                                sliderInput("EWT_l2", HTML("Water content (g cm<sup>-2</sup>)"), min = 0, max = 0.05, value = c(0.01,0.02)),
                                checkboxInput("Dist_ewt_l2", "Use a gauss distribution", value = F),
                                
                                sliderInput("lign_cell2", "Lignin and cellulose content", min = 10, max = 80, value = c(20,60)),
                                checkboxInput("Dist_lign_cell_l2", "Use a gauss distribution", value = F),
                                
                                sliderInput("Nitrogen2", HTML("Nitrogen content (g m<sup>-2</sup>)"), min = 0.3, max = 2,value = c(0.5,1)),
                                checkboxInput("Dist_Nitrogen_l2", "Use a gauss distribution", value = F),
                                
                                sliderInput("cell_d2", HTML("Cell diameter  (m<sup>-6</sup>)"), min = 20, max = 200, value = c(50,100)),
                                checkboxInput("Dist_cell_d_l2", "Use a gauss distribution", value = F),
                                
                                sliderInput("inter_c2", "Intercellular air space", min = 0.01, max = 0.1, value = c(0.05,0.1)),
                                checkboxInput("Dist_inter_c_l2", "Use a gauss distribution", value = F),
                                
                                sliderInput("baseline_abs2", "baseline", min = 0.0004, max = 0.0006, value = c(0.0004,0.0006)),
                                checkboxInput("Dist_baseline_abs_l2", "Use a gauss distribution", value = F),
                                
                                sliderInput("leaf_thick2", "leaf thickness", min = 1, max = 10, value = c(2,5)),
                                checkboxInput("Dist_leaf_thick_l2", "Use a gauss distribution", value = F),
                                
                                sliderInput("albino_abs2", "Albino absorption", min = 0, max = 4, value = c(0,4)),
                                checkboxInput("Dist_albino_abs_l2", "Use a gauss distribution", value = F),
                                
                              ),
                              h4(strrep("-", 60)),
                              h3('canopy parameters :'),
                              # Leaf Model parameters
                              conditionalPanel(
                                condition = "input.canopy_lut == 'fourSAILH'",
                                sliderInput("LAI2", HTML("leaf area index (m<sup>2</sup>/m<sup>-2</sup>)"), min = 0, max = 10,  value = c(0.5,5)),
                                checkboxInput("Dist_lai2", "Use a gauss distribution", value = F),
                                
                                sliderInput("LIDFa2", "LIDFa (°)", min = 0, max = 90, value = c(30,90)),
                                checkboxInput("Dist_lidfa2", "Use a gauss distribution", value = T),
                                
                                sliderInput("hotspot2", "hotspot", min = 0, max = 1, value = c(0.5,0.9)),
                                checkboxInput("Dist_hotspot2", "Use a gauss distribution", value = F),
                                
                                sliderInput("tts2", "tts (deg)", min = 0, max = 90, value = c(0,15)),
                                checkboxInput("Dist_tts2", "Use a gauss distribution", value = T),
                                
                                sliderInput("tto2", "tto (deg)", min = 0, max = 90, value = c(0,15)),
                                checkboxInput("Dist_tto2", "Use a gauss distribution", value = T),
                                
                                sliderInput("psi2", "psi (deg)", min = 0, max = 180, value = c(25,45)),
                                checkboxInput("Dist_psi2", "Use a gauss distribution", value = T),
                                
                                sliderInput("psoil2", "soil factor", min = 0, max = 1, value = c(0.5,0.8)),
                                checkboxInput("Dist_psoil2", "Use a gauss distribution", value = T),
                                
                                
                              ),
                              conditionalPanel(
                                condition = "input.canopy_lut == 'INFORM'",
                                sliderInput("LAI_2", HTML("leaf area index (m<sup>2</sup>/m<sup>-2</sup>)"), min = 0, max = 10,  value = c(0.5,5)),
                                checkboxInput("Dist_lai_2", "Use a gauss distribution", value = T),
                                
                                
                                sliderInput("LIDFa_2", "LIDFa (deg)", min = 0, max = 90, value = c(30,90)),
                                checkboxInput("Dist_lidfa_2", "Use a gauss distribution", value = T),
                                
                                sliderInput("hotspot_2", "hotspot", min = 0, max = 1,value = c(0.5,1)),
                                checkboxInput("Dist_hotspot_2", "Use a gauss distribution", value = F),
                                
                                sliderInput("tts_2", "tts (deg)", min = 0, max = 90, value = c(0,15)),
                                checkboxInput("Dist_tts_2", "Use a gauss distribution", value = T),
                                
                                sliderInput("tto_2", "tto (deg)", min = 0, max = 90, value = c(0,15)),
                                checkboxInput("Dist_tto_2", "Use a gauss distribution", value = T),
                                
                                sliderInput("psi_2", "psi (deg)", min = 0, max = 180, value = c(15,65)),
                                checkboxInput("Dist_psi_2", "Use a gauss distribution", value = F),
                                
                                sliderInput("LAIu_2", HTML("understory LAI (m<sup>2</sup>/m<sup>-2</sup>)"), min = 0.05, max = 2, value = c(0.5,1.5)),
                                checkboxInput("Dist_laiu_2", "Use a gauss distribution", value = F),
                                
                                sliderInput("cd_2", "Crown diameter (m)", min = 0, max = 10, value = c(2,5)),
                                checkboxInput("Dist_cd_2", "Use a gauss distribution", value = F),
                                
                                sliderInput("sd_2",  HTML("Stem density (ha<sup>-1</sup>)"), min = 0, max = 3000,  value = c(500,600)),
                                checkboxInput("Dist_sd_2", "Use a gauss distribution", value = F),
                                
                                sliderInput("h_2", "tree height (m)", min = 10, max = 40,  value = c(25,30)),
                                checkboxInput("Dist_h_2", "Use a gauss distribution", value = F),
                                
                                #  sliderInput("skyl_2", "skyl (fraction-Fixed to 0.1)", min = 0.001, max = 0.4, value = c(0.001,0.4)),
                                sliderInput("psoil_2", "soil factor", min = 0, max = 1,  value = c(0.3,0.6)),
                                checkboxInput("Dist_psoil_2", "Use a gauss distribution", value = T),
                              ),
                              
                            ),
                            mainPanel(
                              h3('Simulated average reflectance'),
                              plotOutput("plot_lut"),
                              
                              #uiOutput("histogram_ui"),
                              actionButton("buttonLUT2", "Generate the LUT"),
                              downloadButton("downloadData2", "Download"),
                              br(),br(),
                              p(em('Note₁: The "Generate the LUT" button is designed to function correctly
                                             only on the first use. Please ensure to click it only once for accurate results.')),
                              p(em('Note₂: Uniform distribution is used as default. Those unrealist simulations were discarded')),
                              p(em('Note₃: Please use the scroll panel to adjust the input ranges, and then generate the Look-Up Table.')),
                              p(em('Note₄: The plant traits shown below represent their respective distributions based on the input parameters.')),
                              # Provide a clear instruction to the user
                              # Add a grey horizontal line
                              tags$hr(style = "border-color: grey;"),
                              p('Please review your selected inputs by choosing between the Min-Max LUT summary or the distribution of plant traits.'),
                              p(em('Note₅:: This option is only active when the simulations are run for the first time; thereafter, it becomes dynamic.')),
                              # Add selectInput to toggle between views
                              selectInput("view_selector", "",
                                          choices = c("Min-Max Inputs (LUT Summary)" = "LUT",
                                                      "Plant Trait Distribution" = "Trait"),
                                          selected = "LUT"), # Set default to LUT summary
                              
                              # Conditional panel for the LUT parameters summary
                              conditionalPanel(
                                condition = "input.view_selector == 'LUT'", # Show LUT summary by default
                                h3('Min-Max Inputs (LUT Summary)'),
                                p("Biochemical, biophysical, structural, and geometric parameters values used for the selected canopy model:"),
                                verbatimTextOutput("params_output") # LUT parameters summary output
                              ),
                              
                              # Conditional panel for the trait distribution plot
                              conditionalPanel(
                                condition = "input.view_selector == 'Trait'", # Show trait distribution plot when selected
                                h3("Plant Trait Distribution"),
                                plotOutput("trait_distribution_plot") # For trait distribution plot
                              ),
                              
                              
                              
                            )## close the mainPanel
                          ) ## close the sidebarLayout
                          
                 ), # Close the tab panel
)
# Define server logic required to draw a histogram ----
server <- function(input, output,session) {
  
  ## 4.1) Personal LUT Section  ---------------------------------------
  
  
  params_lut <- reactive({
    leaf_ranges <- switch(input$leaf_lut,
                          "PROSPECT-PRO" = data.frame(
                            Inputs = c("N", "Cab", "Car", "Anth", "Cbrown", "EWT", "LMA", "Prot", "CBC", "alpha"),
                            Min = c(min(input$N2), min(input$Cab2), min(input$Car2), min(input$Anth2), min(input$Cbrown2),
                                    min(input$EWT2), 0, min(input$Prot2), min(input$CBC2),40),
                            Max = c(max(input$N2), max(input$Cab2), max(input$Car2), max(input$Anth2), max(input$Cbrown2),
                                    max(input$EWT2), 0, max(input$Prot2), max(input$CBC2),40)
                          ),
                          
                          "PROSPECT-D" = data.frame(
                            Inputs = c("N", "Cab", "Car", "Anth", "Cbrown", "EWT", "LMA", "alpha"),
                            Min = c(min(input$N_d2), min(input$Cab_d2), min(input$Car_d2), min(input$Anth_d2), min(input$Cbrown_d2),
                                    min(input$EWT_d2), min(input$LMA_d2), 40),  # Assuming alpha is constant at 40
                            Max = c(max(input$N_d2), max(input$Cab_d2), max(input$Car_d2), max(input$Anth_d2), max(input$Cbrown_d2),
                                    max(input$EWT_d2), max(input$LMA_d2), 40)  # Assuming alpha is constant at 40
                          ),
                          
                          "Liberty" = data.frame(
                            Inputs = c("cell.d", "inter.c", "baseline.abs", "leaf.thick", "albino.abs", "Cab", "EWT", "lign.cell", "Nitrogen"),
                            Min = c(min(input$cell_d2), min(input$inter_c2), 0.0004, min(input$leaf_thick2),
                                    min(input$albino_abs2), min(input$Cab_l2), min(input$EWT_l2), min(input$lign_cell2), min(input$Nitrogen2)),
                            Max = c(max(input$cell_d2), max(input$inter_c2), 0.0004, max(input$leaf_thick2),
                                    max(input$albino_abs2), max(input$Cab_l2), max(input$EWT_l2), max(input$lign_cell2), max(input$Nitrogen2))
                          ),
                          "FLUSPECT-B-Cx" = data.frame(
                            
                            Inputs = c('N',"fqe", "Cab",  "Car", 'Anth',"Cs",  "Cx", "EWT", "LMA","Prot","CBC","alpha"),
                            Min = c(min(input$N_fp_2),min(input$fqe_fp_2), min(input$Cab_fp_2), min(input$Car_fp_2),min(input$Anth_fp_2),
                                    min(input$Cs_fp_2), min(input$Cx_fp_2), min(input$EWT_fp_2),min(input$LMA_fp_2),
                                    min(input$Prot_fp_2),min(input$CBC_fp_2),40),
                            Max = c(max(input$N_fp_2),max(input$fqe_fp_2), max(input$Cab_fp_2), max(input$Car_fp_2),max(input$Anth_fp_2),
                                    max(input$Cs_fp_2), max(input$Cx_fp_2), max(input$EWT_fp_2), max(input$LMA_fp_2),
                                    max(input$Prot_fp_2), max(input$CBC_fp_2),40)
                          ),
                          "FLUSPECT-B" = data.frame(
                            Inputs = c('N',"fqe", "Cab",  "Car", "Cs", "EWT", "LMA", "Cx","alpha"),
                            Min = c(min(input$N_fd_2),min(input$fqe_fd_2), min(input$Cab_fd_2), min(input$Car_fd_2),
                                    min(input$Cs_fd_2), min(input$EWT_fd_2), min(input$LMA_fd_2),min(input$Cx_fd_2), 40),
                            Max = c(max(input$N_fd_2),max(input$fqe_fd_2), max(input$Cab_fd_2), max(input$Car_fd_2),
                                    max(input$Cs_fd_2), max(input$EWT_fd_2), max(input$LMA_fd_2), max(input$Cx_fd_2), 40)
                          )
    )
    
    
    canopy_ranges <- switch(input$canopy_lut,
                            "fourSAILH" = data.frame(
                              Inputs = c("LAI", "TypeLidf", "LIDFa", "LIDFb", "hspot", "tts", "tto", "psi", "psoil"),
                              Min = c(min(input$LAI2), 2, min(input$LIDFa2), 0, min(input$hotspot2),
                                      min(input$tts2), min(input$tto2), min(input$psi2), min(input$psoil2)),
                              Max = c(max(input$LAI2), 2, max(input$LIDFa2), 0, max(input$hotspot2),
                                      max(input$tts2), max(input$tto2), max(input$psi2), max(input$psoil2))
                            ),
                            
                            "INFORM" = data.frame(
                              Inputs = c("LAI", "TypeLidf", "LIDFa", "LIDFb", "hspot", "tts", "tto", "psi", "psoil",
                                         "LAIu", "cd", "sd", "h", "skyl"),
                              Min = c(min(input$LAI_2), 2, min(input$LIDFa_2), 0, min(input$hotspot_2),
                                      min(input$tts_2), min(input$tto_2), min(input$psi_2), min(input$psoil_2),
                                      min(input$LAIu_2), min(input$cd_2), min(input$sd_2), min(input$h_2), 0.1),
                              Max = c(max(input$LAI_2), 2, max(input$LIDFa_2), 0, max(input$hotspot_2),
                                      max(input$tts_2), max(input$tto_2), max(input$psi_2), max(input$psoil_2),
                                      max(input$LAIu_2), max(input$cd_2), max(input$sd_2), max(input$h_2), 0.1)
                            )
    )
    
    return(list(leaf_ranges, canopy_ranges))
  })
  
  
  
  # Reactive expression for the LUT data ---------------------------------------------------
  lut_datasets<- reactive({
    # Get the selected leaf and canopy parameters from the reactive `params` expression
    leaf_params <- params_lut()[[1]]
    canopy_params <- params_lut()[[2]]
    lut_ranges<-rbind(leaf_params,canopy_params)
    return(lut_ranges)
  })
  
  
  # Render the canopy parameters
  output$params_output <- renderPrint({
    # Call the reactive expression
    lut_to_show <- lut_datasets()
    # Print the canopy parameters
    print(lut_to_show)
  })
  
  
  observe({
    # Print the reactive parameters to the console when they change
    params <- params_lut()
    showNotification("Look-up table configured successfully.", type = "message")
    # Print leaf parameters
    # print("Leaf Parameters:")
    # print(params[[1]])  # First element in the list: leaf_params
    
    # Print canopy parameters
    #print("Canopy Parameters:")
    #print(params[[2]])  # Second element in the list: canopy_params
    # Print canopy parameters
    #print("All Parameters:")
    #print(lut_datasets())  # Second element in the list: canopy_params
    
  })
  
  
  w <- reactiveValues(data = NULL)
  # Isolate simulation logic in an eventReactive expression
 
  observeEvent(input$buttonLUT2, {
    
    showNotification("Generating simulations based on the pre-configured Look-Up Table. Please wait...", type = "warning")
    output$plot_lut <- renderPlot({
      
      # Call the Simulations
      data <- ToolsRTM::dataSpec_PDB
      Rsoil.dry  <- data[,11]  # rsoil1 = dry soil
      Rsoil.wet <- data[,12]  # rsoil2 = wet soil
      
      # Show the modal window
      show_modal_spinner()
      
      db.ranges <-lut_datasets()
      print(db.ranges)
      # Determine the distribution based on the checkbox status
      data_lut<- get.LUTfromRanges(LUT=db.ranges,nLUT = input$n_samples_lut, setseed = input$seed_lut,
                                   leaf.model = input$leaf_lut ,
                                   canopy.model = input$canopy_lut,
                                   distribution = 'uniform')
      
      
      # Create a list to store distribution types for each parameter
      distributions <- list(
        ##PROSPECT-PRO
        Cab = ifelse(input$Dist_cab2, 'gauss', 'uniform'),
        Car = ifelse(input$Dist_car2, 'gauss', 'uniform'),
        Anth = ifelse(input$Dist_anth2, 'gauss', 'uniform'),
        Cbrown = ifelse(input$Dist_cbrown2, 'gauss', 'uniform'),
        N = ifelse(input$Dist_n2, 'gauss', 'uniform'),
        EWT = ifelse(input$Dist_ewt2, 'gauss', 'uniform'),
        Prot = ifelse(input$Dist_prot2, 'gauss', 'uniform'),
        CBC = ifelse(input$Dist_cbc2, 'gauss', 'uniform'),
        ##PROSPECT-D
        Cab = ifelse(input$Dist_cab_d2, 'gauss', 'uniform'),
        Car = ifelse(input$Dist_car_d2, 'gauss', 'uniform'),
        Anth = ifelse(input$Dist_anth_d2, 'gauss', 'uniform'),
        Cbrown = ifelse(input$Dist_cbrown_d2, 'gauss', 'uniform'),
        N = ifelse(input$Dist_n_d2, 'gauss', 'uniform'),
        EWT = ifelse(input$Dist_ewt_d2, 'gauss', 'uniform'),
        LMA = ifelse(input$Dist_lma_d2, 'gauss', 'uniform'),
        
        
        ##FLUSPECT-B
        Cab = ifelse(input$Dist_cab_fd_2, 'gauss', 'uniform'),
        Car = ifelse(input$Dist_car_fd_2, 'gauss', 'uniform'),
        fqe = ifelse(input$Dist_fqe_fd_2, 'gauss', 'uniform'),
        Cs = ifelse(input$Dist_cs_fd_2, 'gauss', 'uniform'),
        Cx = ifelse(input$Dist_cx_fd_2, 'gauss', 'uniform'),
        N = ifelse(input$Dist_n_fd_2, 'gauss', 'uniform'),
        EWT = ifelse(input$Dist_ewt_fd_2, 'gauss', 'uniform'),
        LMA = ifelse(input$Dist_lma_fd_2, 'gauss', 'uniform'),
        
        ##FLUSPECT-B-Cx
        Cab = ifelse(input$Dist_cab_fp_2, 'gauss', 'uniform'),
        Car = ifelse(input$Dist_car_fp_2, 'gauss', 'uniform'),
        Anth = ifelse(input$Dist_ant_fp_2, 'gauss', 'uniform'),
        fqe = ifelse(input$Dist_fqe_fp_2, 'gauss', 'uniform'),
        Cs = ifelse(input$Dist_cs_fp_2, 'gauss', 'uniform'),
        Cx = ifelse(input$Dist_cx_fp_2, 'gauss', 'uniform'),
        N = ifelse(input$Dist_n_fp_2, 'gauss', 'uniform'),
        EWT = ifelse(input$Dist_ewt_fp_2, 'gauss', 'uniform'),
        LMA = ifelse(input$Dist_lma_fp_2, 'gauss', 'uniform'),
        Prot = ifelse(input$Dist_prot_fp_2, 'gauss', 'uniform'),
        CBC = ifelse(input$Dist_cbc_fp_2, 'gauss', 'uniform'),
        
        ##Liberty
        cell.d = ifelse(input$Dist_cell_d_l2, 'gauss', 'uniform'),
        inter.c = ifelse(input$Dist_inter_c_l2, 'gauss', 'uniform'),
        baseline.abs = ifelse(input$Dist_baseline_abs_l2, 'gauss', 'uniform'),
        leaf.thick = ifelse(input$Dist_leaf_thick_l2, 'gauss', 'uniform'),
        albino.abs = ifelse(input$Dist_albino_abs_l2, 'gauss', 'uniform'),
        Cab = ifelse(input$Dist_cab_l2, 'gauss', 'uniform'),
        EWT = ifelse(input$Dist_ewt_l2, 'gauss', 'uniform'),
        lign.cell = ifelse(input$Dist_lign_cell_l2, 'gauss', 'uniform'),
        Nitrogen = ifelse(input$Dist_Nitrogen_l2, 'gauss', 'uniform'),
        
        ##fourSAILH
        LAI = ifelse(input$Dist_lai2, 'gauss', 'uniform'),
        LIDFa = ifelse(input$Dist_lidfa2, 'gauss', 'uniform'),
        hspot = ifelse(input$Dist_hspot2, 'gauss', 'uniform'),
        tts = ifelse(input$Dist_tts2, 'gauss', 'uniform'),
        tto = ifelse(input$Dist_tto2, 'gauss', 'uniform'),
        psi = ifelse(input$Dist_psi2, 'gauss', 'uniform'),
        psoil = ifelse(input$Dist_psoil2, 'gauss', 'uniform'),
        
        ##INFOMR
        LAI = ifelse(input$Dist_lai_2, 'gauss', 'uniform'),
        LIDFa = ifelse(input$Dist_lidfa_2, 'gauss', 'uniform'),
        hspot = ifelse(input$Dist_hotspot_2, 'gauss', 'uniform'),
        tts = ifelse(input$Dist_tts_2, 'gauss', 'uniform'),
        tto = ifelse(input$Dist_tto_2, 'gauss', 'uniform'),
        psi = ifelse(input$Dist_psi_2, 'gauss', 'uniform'),
        psoil = ifelse(input$Dist_psoil_2, 'gauss', 'uniform'),
        LAIu = ifelse(input$Dist_laiu_2, 'gauss', 'uniform'),
        cd = ifelse(input$Dist_cd_2, 'gauss', 'uniform'),
        sd = ifelse(input$Dist_sd_2, 'gauss', 'uniform'),
        h = ifelse(input$Dist_h_2, 'gauss', 'uniform') )
      
      # Define expected inputs for each model
      expected_inputs <- list(
        "PROSPECT-D" = c("N", "Cab", "Car", "Anth", "Cbrown", "EWT", "LMA", "alpha"),
        "PROSPECT-PRO" =  c("N", "Cab", "Car", "Anth", "Cbrown", "EWT", "LMA", "Prot", "CBC", "alpha"),
        'FLUSPECT-B' = c('N',"fqe", "Cab",  "Car", "Cs", "EWT", "LMA", "Cx","alpha"),
        'FLUSPECT-B-Cx' = c('N',"fqe", "Cab",  "Car", 'Anth',"Cs",  "Cx", "EWT", "LMA","Prot","CBC","alpha"),
        "Liberty" =   c("cell.d", "inter.c", "baseline.abs", "leaf.thick", "albino.abs", "Cab", "EWT", "lign.cell", "Nitrogen"),
        "fourSAILH" = c("LAI", "TypeLidf", "LIDFa", "LIDFb", "hspot", "tts", "tto", "psi", "psoil"),
        "INFORM" = c("LAI", "TypeLidf", "LIDFa", "LIDFb", "hspot", "tts", "tto", "psi", "psoil",
                     "LAIu", "cd", "sd", "h", "skyl") )
      
      # Check if the selected leaf or canopy model requires special handling
      if (input$leaf_lut %in% names(expected_inputs) || input$canopy_lut %in% names(expected_inputs)) {
        # Get expected parameters for the selected model
        model_inputs <- expected_inputs[[input$leaf_lut]] %||% expected_inputs[[input$canopy_lut]]
        
        # Apply distributions only for the specified parameters present in the data
        for (param in model_inputs) {
          if (param %in% names(distributions) && distributions[[param]] == 'gauss') {
            # Calculate standard deviation for Gaussian distribution
            sd_ <- (max(data_lut[[param]]) - min(data_lut[[param]])) / 6  # Assuming 99.7% coverage for min/max in normal distribution
            data_lut[[param]] <- stats::rnorm(input$n_samples_lut, mean = mean(data_lut[[param]], na.rm = TRUE), sd = sd_)
          }
        }
      }
      # We consider only LIDFa
      data_lut$TypeLidf <- 2
      #Adding alpha to 40 (default in PROSPECT)
      data_lut$alpha <- 40
      # Alternatively, if you want to return it as a data frame in a reactive value
      w$data <- data_lut  # Store LUT in reactive values for later use
      
      
      if ((input$leaf_lut == 'PROSPECT-PRO') & (input$canopy_lut == 'fourSAILH')) {
        
        # choose number of processors/cores
        no_cores <- parallel::detectCores() - 2
        cl <- parallel::makeCluster(no_cores)
        doParallel::registerDoParallel(cl)
        sim.rfl<-list()
        
        sims<-foreach::foreach(i=1:input$n_samples_lut) %dopar% {
        
          
          psoil = data_lut[i,'psoil']
          rsoil_<- c(psoil*Rsoil.dry+(1-psoil)*Rsoil.wet)
          
          data.foursail<-  ToolsRTM::foursail(inputLUT=data_lut[i,],rsoil=rsoil_,LeafModel = 'PROSPECT-PRO')
          rdot<-data.foursail[[1]]
          rsot<-data.foursail[[2]]
          rfl.prosail<- ToolsRTM::Compute_BRF(rdot=rdot,rsot=rsot,tts=data_lut[i,'tts'],data.light= ToolsRTM::dataSpec_PDB)
          sim.rfl[[i]]<-rfl.prosail
          
        } ##end paralle
        
        # Close the cluster after parallel processing is complete
        parallel::stopCluster(cl)
        showNotification("LUT with simulations done successfully.", type = "message")
        
        
      } else if ((input$leaf_lut == 'PROSPECT-D') & (input$canopy_lut == 'fourSAILH')) {
        
        ## choose number of processors/cores
        no_cores <- parallel::detectCores() - 2
        cl <- parallel::makeCluster(no_cores)
        doParallel::registerDoParallel(cl)
        sim.rfl<-list()
        sims<-foreach::foreach(i=1:input$n_samples_lut) %dopar% {
        
          psoil = data_lut[i,'psoil']
          rsoil_<- c(psoil*Rsoil.dry+(1-psoil)*Rsoil.wet)
          
          data.foursail<- ToolsRTM::foursail(inputLUT=data_lut[i,],rsoil=rsoil_,LeafModel = 'PROSPECT-D')
          rdot<-data.foursail[[1]]
          rsot<-data.foursail[[2]]
          rfl.prosail<- ToolsRTM::Compute_BRF(rdot=rdot,rsot=rsot,tts=data_lut[i,'tts'],data.light= ToolsRTM::dataSpec_PDB)
          sim.rfl[[i]]<-rfl.prosail
          
        }
        
        # Close the cluster after parallel processing is complete
        parallel::stopCluster(cl)
        showNotification("LUT with simulations done successfully.", type = "message")
        
      } else if ((input$leaf_lut == 'FLUSPECT-B') & (input$canopy_lut == 'fourSAILH')) {
        
        ## choose number of processors/cores
        no_cores <- parallel::detectCores() - 2
        cl <- parallel::makeCluster(no_cores)
        doParallel::registerDoParallel(cl)
        sim.rfl<-list()
        sims<-foreach::foreach(i=1:input$n_samples_lut) %dopar% {
        
          psoil = data_lut[i,'psoil']
          rsoil_<- c(psoil*Rsoil.dry+(1-psoil)*Rsoil.wet)
          
          data.foursail<- ToolsRTM::foursail(inputLUT=data_lut[i,],rsoil=rsoil_[1:2001],LeafModel = 'Fluspect-B')
          rdot<-data.foursail[[1]]
          rsot<-data.foursail[[2]]
          rfl.prosail<-  ToolsRTM::Compute_BRF(rdot=rdot,rsot=rsot,tts=data_lut[i,'tts'],data.light= ToolsRTM::dataSpec_PDB,short.waves = T)
          sim.rfl[[i]]<-rfl.prosail
          
        }
        
        # Close the cluster after parallel processing is complete
        parallel::stopCluster(cl)
        showNotification("LUT with simulations done successfully.", type = "message")
        
      } else if ((input$leaf_lut == 'FLUSPECT-B-Cx') & (input$canopy_lut == 'fourSAILH')) {
        
        ## choose number of processors/cores
        no_cores <- parallel::detectCores() - 2
        cl <- parallel::makeCluster(no_cores)
        doParallel::registerDoParallel(cl)
        sim.rfl<-list()
        sims<-foreach::foreach(i=1:input$n_samples_lut, .export = c("loadFunctions",'loadRDa')) %dopar% {
          loadRDa("www/data/ToolsRTM")
          loadRDa("www/data")
          loadFunctions()
          psoil = data_lut[i,'psoil']
          rsoil_<- c(psoil*Rsoil.dry+(1-psoil)*Rsoil.wet)
          
          data.foursail<-  ToolsRTM::foursail(inputLUT=data_lut[i,],rsoil=rsoil_[1:2001],LeafModel = 'Fluspect-B-Cx')
          rdot<-data.foursail[[1]]
          rsot<-data.foursail[[2]]
          rfl.prosail<- ToolsRTM::Compute_BRF(rdot=rdot,rsot=rsot,tts=data_lut[i,'tts'],data.light= ToolsRTM::dataSpec_PDB, short.waves = T)
          sim.rfl[[i]]<-rfl.prosail
          
        }
        
        # Close the cluster after parallel processing is complete
        parallel::stopCluster(cl)
        showNotification("LUT with simulations done successfully.", type = "message")
        
      } else if ((input$leaf_lut == 'Liberty') & (input$canopy_lut == 'fourSAILH')) {
        
        ## choose number of processors/cores
        no_cores <- parallel::detectCores() - 2
        cl <- parallel::makeCluster(no_cores)
        doParallel::registerDoParallel(cl)
        sim.rfl<-list()
        sims<-foreach::foreach(i=1:input$n_samples_lut) %dopar% {
          
          psoil = data_lut[i,'psoil']
          rsoil_<- c(psoil*Rsoil.dry+(1-psoil)*Rsoil.wet)
          
          data.foursail<-  ToolsRTM::foursail(inputLUT=data_lut[i,],rsoil=rsoil_,LeafModel = 'Liberty')
          rdot<-data.foursail[[1]]
          rsot<-data.foursail[[2]]
          rfl.prosail<- ToolsRTM::Compute_BRF(rdot=rdot,rsot=rsot,tts=data_lut[i,'tts'],data.light= ToolsRTM::dataSpec_PDB)
          sim.rfl[[i]]<-rfl.prosail
          
        } ##end paralle
        
        # Close the cluster after parallel processing is complete
        parallel::stopCluster(cl)
        showNotification("LUT with simulations done successfully.", type = "message")
        
        
      } else if ((input$leaf_lut == 'PROSPECT-PRO') & (input$canopy_lut == 'INFORM')) {
        
        ## choose number of processors/cores
        no_cores <- parallel::detectCores() - 2
        cl <- parallel::makeCluster(no_cores)
        doParallel::registerDoParallel(cl)
        sim.rfl<-list()
        sims<-foreach::foreach(i=1:input$n_samples_lut) %dopar% {
    
          
          psoil = data_lut[i,'psoil']
          rsoil_<- c(psoil*Rsoil.dry+(1-psoil)*Rsoil.wet)
          
          rfl.inform<-  ToolsRTM::inform(inputLUT=data_lut[i,],rsoil=rsoil_,LeafModel = 'PROSPECT-PRO')
          sim.rfl[[i]]<-rfl.inform
          
        } ##end paralle
        
        # Close the cluster after parallel processing is complete
        parallel::stopCluster(cl)
        showNotification("LUT with simulations done successfully.", type = "message")
        
        
      } else if ((input$leaf_lut == 'PROSPECT-D') & (input$canopy_lut == 'INFORM')) {
        
        ## choose number of processors/cores
        no_cores <- parallel::detectCores() - 2
        cl <- parallel::makeCluster(no_cores)
        doParallel::registerDoParallel(cl)
        sim.rfl<-list()
        sims<-foreach::foreach(i=1:input$n_samples_lut) %dopar% {
       
          psoil = data_lut[i,'psoil']
          rsoil_<- c(psoil*Rsoil.dry+(1-psoil)*Rsoil.wet)
          
          rfl.inform<- ToolsRTM::inform(inputLUT=data_lut[i,],rsoil=rsoil_,LeafModel = 'PROSPECT-D')
          sim.rfl[[i]]<-rfl.inform
          
        } ##end paralle
        # Close the cluster after parallel processing is complete
        parallel::stopCluster(cl)
        showNotification("LUT with simulations done successfully.", type = "message")
        
      } else if ((input$leaf_lut == 'FLUSPECT-B') & (input$canopy_lut == 'INFORM')) {
        
        ## choose number of processors/cores
        no_cores <- parallel::detectCores() - 2
        cl <- parallel::makeCluster(no_cores)
        doParallel::registerDoParallel(cl)
        sim.rfl<-list()
        sims<-foreach::foreach(i=1:input$n_samples_lut) %dopar% {
        
          psoil = data_lut[i,'psoil']
          rsoil_<- c(psoil*Rsoil.dry+(1-psoil)*Rsoil.wet)
          
          rfl.inform<- ToolsRTM::inform(inputLUT=data_lut[i,],rsoil=rsoil_[1:2001],LeafModel = 'Fluspect-B')
          sim.rfl[[i]]<-rfl.inform
          
        } ##end paralle
        # Close the cluster after parallel processing is complete
        parallel::stopCluster(cl)
        showNotification("LUT with simulations done successfully.", type = "message")
        
      } else if ((input$leaf_lut == 'FLUSPECT-B-Cx') & (input$canopy_lut == 'INFORM')) {
        
        ## choose number of processors/cores
        no_cores <- parallel::detectCores() - 2
        cl <- parallel::makeCluster(no_cores)
        doParallel::registerDoParallel(cl)
        sim.rfl<-list()
        sims<-foreach::foreach(i=1:input$n_samples_lut) %dopar% {
        
          psoil = data_lut[i,'psoil']
          rsoil_<- c(psoil*Rsoil.dry+(1-psoil)*Rsoil.wet)
          
          rfl.inform<- ToolsRTM::inform(inputLUT=data_lut[i,],rsoil=rsoil_[1:2001],LeafModel = 'Fluspect-B-Cx')
          sim.rfl[[i]]<-rfl.inform
          
        } ##end paralle
        # Close the cluster after parallel processing is complete
        parallel::stopCluster(cl)
        showNotification("LUT with simulations done successfully.", type = "message")
        
      } else if ((input$leaf_lut == 'Liberty') & (input$canopy_lut == 'INFORM')) {
        
        
        ## choose number of processors/cores
        no_cores <- parallel::detectCores() - 2
        cl <- parallel::makeCluster(no_cores)
        doParallel::registerDoParallel(cl)
        sim.rfl<-list()
        sims<-foreach::foreach(i=1:input$n_samples_lut) %dopar% {
    
          psoil = data_lut[i,'psoil']
          rsoil_<- c(psoil*Rsoil.dry+(1-psoil)*Rsoil.wet)
          rfl.inform<-  ToolsRTM::inform(inputLUT=data_lut[i,],rsoil=rsoil_,LeafModel = 'Liberty')
          sim.rfl[[i]]<-rfl.inform
          
        } ##end paralle
        
        # Close the cluster after parallel processing is complete
        parallel::stopCluster(cl)
        showNotification("LUT with simulations done successfully.", type = "message")
        
      }
      
      
      sim.canopy<-do.call(rbind,sims)
      if ((input$leaf_lut == 'FLUSPECT-B') | (input$leaf_lut == 'FLUSPECT-B-Cx') ){
        wave<-seq(400,2500,1)[1:2001]
      } else{
        wave<-seq(400,2500,1)
      }
      
      
      IDs<-c(1:input$n_samples_lut)
      # Convert matrix to data frame
      df <- data.frame(sim.canopy)
      print(dim(df))
      df$row <- 1:nrow(df)  # Add a row identifier
      #df <-df[!complete.cases(df), ]
      df <-na.omit(df)
      # Extract the row indices that were kept in df_clean
      rows_kept <- df$row
      
      
      print('removing cases')
      print(dim(df))
      print(any(is.na(df)))
      #print(colSums(is.na(df)))
      
      if (any(is.na(df)) == T){
        
        print(dim( df[!complete.cases(df), ]))
      }
      
      ###
      if (input$sensor_lut == 'RTM') {
        
        df <-na.omit(df)
        # Now, select the same rows from data_lut based on the remaining rows in df
        data_lut_filtered <- data_lut[rownames(df), ]
        IDs<-c(1:input$n_samples_lut)[rownames(df)]
        print(get.plots(df=df,wave=wave))
        df.to_export <- cbind(IDs,data_lut_filtered,df)
        
      } else {
        
        df.sensor <- subset(sensor.characteristics, Sensor == input$sensor_lut)
        color.selected <- c('forestgreen','red','darkgoldenrod2','dodgerblue3','darkseagreen3','indianred2')
        print(df.sensor)
        
        
        if (input$sensor_lut == 'Sentinel2a'){
          color.i <-color.selected[1]
          Band_ <- paste('B',df.sensor$channel,sep='')
          Band_ <- c('B1','B2','B3','B4','B5','B6','B7','B8','B8A','B9','B10','B11','B12')
        } else if (input$sensor_lut == 'Sentinel2b'){
          color.i <-color.selected[2]
          Band_ <- paste('B',df.sensor$channel,sep='')
          Band_ <- c('B1','B2','B3','B4','B5','B6','B7','B8','B8A','B9','B10','B11','B12')
        } else if (input$sensor_lut == 'Landsat4'){
          color.i <-color.selected[3]
          Band_ <- paste('B',df.sensor$channel,sep='')
          Band_ <- c('B1','B2','B3','B4','B5','B7')
        } else if (input$sensor_lut == 'Landsat5'){
          color.i <-color.selected[4]
          Band_ <- paste('B',df.sensor$channel,sep='')
          Band_ <- c('B1','B2','B3','B4','B5','B7')
        } else if (input$sensor_lut == 'Landsat7'){
          color.i <-color.selected[5]
          Band_ <- paste('B',df.sensor$channel,sep='')
          Band_ <- c('B1','B2','B3','B4','B5','B7')
        } else if (input$sensor_lut == 'Landsat8'){
          color.i <-color.selected[6]
          Band_ <- paste('B',df.sensor$channel,sep='')
          Band_ <- c('B1','B2','B3','B4','B9','B9','B6','B7')
        }
        
        fwhm <-(df.sensor$lb -df.sensor$ub)
        center_wvl <-c(df.sensor$average)
        
        
        # Use these indices to remove the same rows from sim.canopy
        sim.canopy <- sim.canopy[rows_kept, ]
        
        # Initialize an empty list to store interpolated reflectance data
        Sz_ <-list()
        Sz_df <- list()
        # Perform linear interpolation for each row in sim.canopy
        for (i in 1:nrow(sim.canopy)) {
          # Perform linear interpolation to resample to Sentinel-2 bands
          reflectance.i <- signal::interp1(wave, sim.canopy[i, ], center_wvl, method = "spline")
          # Store the interpolated reflectance in the list
          Sz_df[[i]] <- data.frame(center_wvl = center_wvl, reflectance = reflectance.i)
          Sz_[[i]] <-reflectance.i
        }
        
        # Convert the list to a matrix
        df.Sz <- do.call(rbind, Sz_df)
        df.Sa_rfl <- as.data.frame(do.call(rbind, Sz_))
        colnames(df.Sa_rfl) <-Band_
        IDs <-rownames(df.Sa_rfl)  # Change this to ensure the IDs match df.Sa_rfl's dimensions
        
        df.to_export <- cbind(IDs,data_lut[rows_kept,],df.Sa_rfl)
        
        # Calculate average, 25th percentile, and 50th percentile for each band
        stats.SE2 <- df.Sz %>%
          group_by(center_wvl) %>%
          summarise(
            average = mean(reflectance,na.rm=T),
            median = median(reflectance,na.rm=T),
            percentile_25 = quantile(reflectance, 0.25,na.rm=T),
            percentile_50 = quantile(reflectance, 0.50,na.rm=T),
            percentile_75 = quantile(reflectance, 0.75,na.rm=T))
        # Plot using ggplot2
        plot_ <-ggplot(stats.SE2, aes(x = center_wvl)) + #ylim(0,0.8) +
          xlim(400,2500) +
          geom_line(aes(y = average), color = color.i, size = 0.8) +
          geom_line(aes(y = median),  linetype = "dashed", color = color.i, size = 0.8) +
          geom_ribbon(aes(ymin = percentile_25, ymax = percentile_75), linetype = "dashed",fill = color.i, alpha = 0.3) +
          labs(
            title = "",
            x = "wavelength (nm)",
            y = "Reflectance"
          ) +
          theme_bw() +  theme(
            text = element_text(size = 14, face='bold'),  # Increase the text size
            axis.title = element_text(size = 16, face = "bold"),  # Make axis titles bold
            axis.text = element_text(face = "bold"),  # Make axis numbers bold
            plot.title = element_text(face = "bold"),  # Make plot title bold
            plot.subtitle = element_text(size = 14, face='bold')  # Adjust subtitle size
          )
        showNotification("LUT successfully adjusted to the spectral resolution of the selected satellite.", type = "message")
        
        print(plot_)
        
      }
      remove_modal_spinner()
      
      #Create a reactive data frame for saving data
      save_data <- reactive({
        df.to_export
      })
      # Add the download handler for saving data
      output$downloadData2 <- downloadHandler(
        filename = function() {
          paste("LUT_", input$leaf_lut,'_', input$canopy_lut,'_', input$sensor_lut, ".csv", sep = "")
        },
        content = function(file) {
          write.csv(save_data(), file, row.names = FALSE)
        }
      )
      
    }) #end outplots
    
  })    #end observer pannel
  
  
  # Print LUT data outside the observeEvent
  observe({
    req(w$data)  # Ensure w$data is not NULL
    print(head(w$data))  # Print the LUT data frame to the console
  })
  
  output$trait_distribution_plot <- renderPlot({
    req(w$data)
    # Filter out traits where all values are the same
    filtered_data <- w$data %>%
      dplyr::select(where(~ sd(.) != 0))  # Keep only columns with non-zero standard deviation
    
    # Convert data to long format
    data_long <- filtered_data %>%
      pivot_longer(cols = everything(), names_to = "Trait", values_to = "Value")
    
    # Create the ggplot
    ggplot(data_long, aes(x = Value)) +
      geom_histogram(bins = 30, fill = "steelblue", color = "black") +
      facet_wrap(~Trait, scales = "free") +
      theme_bw() +
      labs(title = "",x = "", y = "Frequency") + theme(
        text = element_text(size = 14, face='bold'),  # Increase the text size
        axis.title = element_text(size = 16, face = "bold"),  # Make axis titles bold
        axis.text = element_text(face = "bold"),  # Make axis numbers bold
        plot.title = element_text(face = "bold"),  # Make plot title bold
        plot.subtitle = element_text(size = 14, face='bold')  # Adjust subtitle size
      )
  })
  
  # Alternatively,  show the first few rows of the data
  output$lut_head <- renderTable({
    req(w$data)  # Ensure w$data is not NULL
    head(w$data)  # Show the first few rows of the data frame
  })
  
  
}

# Create Shiny object
shinyApp(ui = ui, server = server)

