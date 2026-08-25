

rm(list= ls())

# 1. load the main libraries  -----

library(ToolsRTM)
required_packages <- c("shiny", "shinythemes", "ggplot2", "dplyr",'DT')

# Check for missing packages and install them if necessary
missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}
# Load the libraries
lapply(required_packages, library, character.only = TRUE)

# Define UI
ui <- navbarPage("Soil Reflectance Simulator", theme = shinytheme("flatly"),

  # Main Tab Panel for Running the Simulation
  tabPanel("MARMIT Model",
           sidebarLayout(
             sidebarPanel(
               h4("MARMIT version"),
               selectInput("version_marmit", "Version", choices = list(
                 "MARMIT 1" = "marmit1",

                 "MARMIT 2" = "marmit2"
               )),

               h4("Soil database"),
               selectInput("database", "Database", choices = list(
                 "Bablet 2016" = "Bablet_2016",
                 "Humper 2015" = "Humper_2015",
                 "Dupiau-2020" = "Dupiau_2020",
                 "Lesaignoux 2008" = "Lesaignoux_2008",
                 "Lobell 2002" = "Lobell_2002",
                 "Marcq 2012" = "Marcq_2012",
                 "Liu 2002" = "Liu_2002",
                 "Philpot 2014" = "Philpot_2014"
               )),
               sliderInput("id", "Soil ID", min = 1, max = 17, step = 1, value = 1),


               # Conditional panel for MARMIT 1
               conditionalPanel(
                 condition = "input.version_marmit == 'marmit1'",
                 h4("Soil Parameters"),
                 sliderInput("eps_marmit1", "Wet soil surface ratio (ε)", min = 0, max = 1, step = 0.05,value = 0.3),
                 sliderInput("L_marmit1", "Thickness of water layer (cm)", min = 0.0001, max = 0.5, step = 0.001, value = 0.05)
               ),
               # Conditional panel for MARMIT 2
               conditionalPanel(
                 condition = "input.version_marmit == 'marmit2'",
                 h4("Soil Parameters"),
                 sliderInput("eps_marmit2", "Wet soil surface ratio (ε)", min = 0, max = 1,  step = 0.05, value = 0.3),
                 sliderInput("L_marmit2", "Thickness of water layer (cm)", min = 0.0001, max = 0.5, step = 0.001, value = 0.05),
                 sliderInput("d_i", "Volume fraction of soil particles", min = 0.0001, max = 0.1, step = 0.005, value = 0.0005)
               ),
               checkboxInput(inputId = "include_soil_reference",
                             label = "Include soil dry reference",
                             value = F) , # Default to TRUE or FALSE as per your preference
               downloadButton("downloadData", "Download MARMIT Data")

             ),
             mainPanel(
               tabsetPanel(
                 id = "main_tabs",
                 tabPanel("Soil Simulation ",
                          h4(""),
                          plotOutput("plot_marmit", height = "400px"),
                          DT::DTOutput("table_marmit")
                 ),

                 tabPanel("Model Information",
                          h4("MARMIT Model Overview"),
                          p("MARMIT (Multilayer rAdiative tRansfer Model of soIl reflecTance) is a radiative transfer model
             that predicts the spectral reflectance of bare soil from 400 nm to 2500 nm with a 1 nm step in the solar domain.
             The model estimates reflectance based on surface water content and various soil parameters."),
                          p("MARMIT-2 is the improved version, providing more accurate predictions of soil reflectance."),

                          h4("References"),
                          p("[1] Bablet A., Vu P.V.H., Jacquemoud S., Viallefont-Robinet F., Fabre S., Briottet X., Sadeghi M., Whiting M.L.,
             Baret F., and Tian J. (2018), MARMIT: a multilayer radiative transfer model of soil reflectance to estimate surface soil moisture
             content in the solar domain (400–2500 nm), Remote Sensing of Environment, 217:1-17.
             ", a("https://doi.org/10.1016/j.rse.2018.07.031", href="https://doi.org/10.1016/j.rse.2018.07.031")),
                          p("[2] Dupiau A., Jacquemoud S., Briottet X., Fabre S., Viallefont-Robinet F., Philpot W., Di Biagio C., Hébert H., and
             Formenti P. (2022), MARMIT-2: an improved version of the MARMIT model to predict soil reflectance as a function of surface
             water content in the solar domain, Remote Sensing of Environment, 272:112951.
             ", a("https://doi.org/10.1016/j.rse.2022.112951", href="https://doi.org/10.1016/j.rse.2022.112951")),
                          p("More information and the GitLab repository for the MARMIT model can be found ",
                            a("here", href="https://pss-gitlab.math.univ-paris-diderot.fr/marmit/marmit"), ".")

                 ),
               )
             )
           )
  ),


)
# Water-optics reference tables never depend on any input (database, soil ID,
# L, eps, ...) -- load them once at app startup instead of inside a reactive
# that re-reads them from disk on every slider drag (L/eps sliders trigger
# ~10x/second while dragging; these are the two largest bundled CSVs).
n_w_full <- read.csv("www/marmit/parameters/n_segelstein.csv", sep = "\t")
alpha_w_full <- read.csv("www/marmit/parameters/alpha_buikouwie.csv", sep = "\t")

# Define server logic
server <- function(input, output, session) {

  # Dynamically update the 'id' slider based on the selected database
  observeEvent(input$database, {
    db_limits <- list(
      "Bablet_2016" = 17, "Dupiau_2020" = 8, "Humper_2015" = 57,
      "Lesaignoux_2008" = 32, "Liu_2002" = 92, "Lobell_2002" = 4,
      "Marcq_2012" = 9, "Philpot_2014" = 3
    )
    max_id <- db_limits[[input$database]]
    updateSliderInput(session, "id", min = 1, max = max_id)
  })

  # Reactive to store simulation results
  sim_results <- reactiveValues(spectral = NULL, params = NULL, dry_soil =NULL)

  # Dynamically update d_i based on the selected MARMIT version
  d_i <- reactive({
    if (input$version_marmit == "marmit1") {
      NULL  # d_i is NULL for MARMIT 1
    } else {
      input$d_i  # d_i is taken from input for MARMIT 2
    }
  })

  # Everything that only depends on the chosen database/soil ID (not on the
  # L/eps sliders) -- separated into its own reactive so it's cached and only
  # re-runs when input$database/input$id actually change, not on every slider
  # drag (Shiny's reactive graph tracks exactly which inputs each reactive()
  # reads and only invalidates it when those specific inputs change).
  dry_soil_data <- reactive({
    database <- input$database
    req(database)
    df <- read.csv(paste0("www/marmit/databases/", database, "/", database, ".csv"))

    # Set wavelengths based on the first reflectance spectrum
    R <- read.csv(paste0("www/marmit/databases/", database, "/spectra/", df$Refl_file[1]), sep = "\t")

    wlmin <- max(min(R$Wvl), 400)
    wlmax <- max(R$Wvl)
    wls <- seq(wlmin, wlmax, by = 1)

    # n_w_full/alpha_w_full are loaded once at app startup (see top of file) --
    # just subset the wavelength range here, no disk read.
    n_w <- n_w_full[n_w_full$Wvl >= wlmin & n_w_full$Wvl <= wlmax, "n"]
    alpha_w <- alpha_w_full[alpha_w_full$Wvl >= wlmin & alpha_w_full$Wvl <= wlmax, "alpha"]

    # Read dry soil reflectance
    ID <- input$id
    req(ID)
    df1 <- df[df$ID == ID, ]
    req(nrow(df1) > 0)

    Rd <- read.csv(paste0("www/marmit/databases/", database, "/spectra/", df1$Refl_file[df1$SMCg == min(df1$SMCg)]), sep = "\t")
    Rd <- Rd[Rd$Wvl >= min(wls) & Rd$Wvl <= max(wls), "R"]

    # Fetch sigmoid parameters
    K <- as.numeric(unique(df1$K))
    a <- as.numeric(unique(df1$a))
    psi <- as.numeric(unique(df1$psi))

    list(wls = wls, n_w = n_w, alpha_w = alpha_w, Rd = Rd, K = K, a = a, psi = psi)
  })

  # Runs the (cheap, vector-math-only) MARMIT model itself -- reacts to the
  # L/eps sliders, reusing dry_soil_data()'s cached, already-loaded tables
  # instead of re-reading them from disk.
  observe({
    dsd <- dry_soil_data()
    req(dsd)

    # Soil parameters
    if (input$version_marmit == 'marmit2') {
      L <- input$L_marmit2
      eps <- input$eps_marmit2
    } else {
      L <- input$L_marmit1
      eps <- input$eps_marmit1
    }

    di_value <- d_i()  # Use the reactive value of d_i
    n_i <- 1.53
    k_i <- 0.001

    # Ensure version_marmit has a value before proceeding
    req(input$version_marmit)

    # Run the MARMIT model based on the selected version
    if (input$version_marmit == 'marmit2') {
      Rw <- get.marmit2(dsd$n_w, dsd$alpha_w, n_i, k_i, dsd$Rd, L, eps, di_value, dsd$wls)
    } else {
      Rw <- get.marmit1(dsd$n_w, dsd$alpha_w, dsd$Rd, L, eps)
    }
    # Calculate soil moisture content
    phi <- L * eps
    SMC <- sigmoid.soil(phi, dsd$K, dsd$a, dsd$psi)

    # Store simulation results
    sim_results$spectral <- data.frame(Wavelength = dsd$wls, Reflectance = Rw)
    sim_results$params <- data.frame(L = L, phi = phi, eps = eps, SMCg = SMC)
    sim_results$dry_soil <- data.frame(Wavelength = dsd$wls, Reflectance = dsd$Rd)

    showNotification("soil reflectance done successfully.", type = "message")
  })


  # Create a reactive data frame for saving data
  save_data <- reactive({
    req(sim_results$spectral)
    df.to_export <- sim_results$spectral

    return(df.to_export)  # Ensure to return the dataframe
  })

  # Plot spectral results
  output$plot_marmit <- renderPlot({

    req(sim_results$spectral)

    plot_soil <-ggplot(sim_results$spectral, aes(x = Wavelength, y = Reflectance)) +
      geom_line(linewidth = 1, color = "goldenrod4") +
      labs(x = "Wavelength (nm)", y = "Soil Reflectance") +
      labs(color = "Legend") + theme_bw() +
      theme(legend.position = 'none',
            legend.box.background = element_rect(color = "black",linewidth=1),
            plot.title = element_text(hjust = 0.5, size=14,face="bold"),
            axis.title = element_text(face="bold", size=14),
            legend.text = element_text(face="bold", size=12),
            axis.text.y=element_text(hjust = 0.5, size=12,face="bold"),
            axis.text.x=element_text(hjust = 0.5, size=12,face="bold"),
            legend.title=element_blank())

    # Add the soil dry reflectance in black if the checkbox is checked
    if (input$include_soil_reference) {
      req(sim_results$dry_soil)
      plot_soil <-plot_soil + geom_line(data=sim_results$dry_soil,aes(x = Wavelength, y = Reflectance),
                                        color = "black", linetype = "solid", linewidth = 1)

    }

    # Plot the combined results
    print(plot_soil)

  })

  # Reactive value to track download success
  download_success <- reactiveVal(FALSE)

    # Create a download link
    output$downloadData <- downloadHandler(
      filename = function() {
        paste("soil_spectrum_", Sys.Date(), ".csv", sep = "")  # Added date to filename for uniqueness
      },
      content = function(file) {
        # Get the selected folder path

        write.csv(save_data(), file, row.names = FALSE)
      }
    )
    # Set download success to TRUE
    download_success(TRUE)


  # Observe the download success reactive value to show notification
  observeEvent(download_success(), {
    if (download_success()) {
      showNotification("Downloaded MARMIT data", type = "message")
      # Reset the reactive value to FALSE after showing the notification
      download_success(FALSE)
    }
  })


  # Show input parameters in a DT table with only 1 row and no search functionality
  output$table_marmit <- DT::renderDT({
    req(sim_results$params)
    datatable(sim_results$params[1, ], options = list(
      dom = 't',  # Only show the table (no search, no pagination, no info)
      paging = FALSE  # Disable pagination
    ))
  })


}
# Create the Shiny app object
shinyApp(ui = ui, server = server)
