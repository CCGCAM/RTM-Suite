
rm(list= ls())

# 1. load the main libraries  -----

library(ToolsRTM)
required_packages <- c("shiny", 'shinybusy',"shinythemes", "ggplot2","reshape2","htmlwidgets",
                       'foreach','doParallel','parallel', "dplyr",'DT','tidyr','shinyWidgets',
                       'sf','leaflet','fs','rstac','gdalcubes','terra','leafem')

# 2.Check for missing packages and install them if necessary  -----

missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}
# Load the libraries
lapply(required_packages, library, character.only = TRUE)

# 3.Define UI -----

ui <- navbarPage("Sentinel-2's Scenario ",theme = shinytheme("flatly"),

                 tabPanel(title = "STAC application",
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    # Sidebar panel for inputs ----

    # Sidebar panel for inputs ----
    sidebarPanel(width = 3,

                 selectInput("scenario", "Select Scenario",
                             choices = c("Bark beetle outbreak", "San Rossore")),  # Scenarios choices
                 selectInput("server_map", "Select Spatial platform",
                             choices = c("Microsoft Planetary Computer",
                                     "Amazon Web Service")),  # Platforms for spatial data access

                 selectInput("satellite_collection", "Satellite Collections:",
                             choices = c("Sentinel-2 Collection" )),
                 selectInput("aggregation_method", "Aggregation methods:",
                             choices = c("Minimum","Maximum", "Mean", "Median" ), #'First'
                             selected = "Median"),


                 selectInput("resampling_method", "Resampling methods:",
                             choices = c("Near","Bilinear", "Bicubic" ),
                             selected = "Bilinear"),
                 numericInput("bbox_cube", "Composite window (m)", min = 20, max = 10000, value = 600, step = 10),

            #     dateRangeInput("date_range_collection", "Select a time period:",
             #                   start = "2023-06-01", end = "2023-06-15",
              #                  min = "2020-06-23", max = "2024-06-23"),

                 sliderInput("cloud_collection", "Maximum Cloud Coverage (%):", min = 5, max = 100, value = 25, step = 1),


                 fluidRow(
                   column(6,  # First column for the action button
                          actionButton("get_cube", "Get the collection")
                   ),

                   column(6,  # Second column for the conditional download button
                          conditionalPanel(
                            condition = "output.raster_ready == true",  # Only show if the output is ready
                            downloadButton("download_tiff", label = "Download")
                          )
                   )
                 )
    ),

    # Main panel for displaying outputs ----
    mainPanel(width = 6,

              # Create a tabset panel
              tabsetPanel(id = "tabs",
                tabPanel("Map viewer",
                         br(),
                         p(style = "text-align: justify;", HTML('This app allows you to retrieve satellite imagery via the Microsoft Planetary STAC API to analyze vegetation health over time. Start by selecting a scenario, time period, and a buffer area around your chosen polygon’s centroid.
                         Then, click on the map to view
                         spectral profiles for the selected image.')),
                         # Add custom CSS styling for the map border and slider
                         tags$head(
                           tags$style(HTML("
                           #map {
                           border: 2px solid #007BFF; /* Blue border */
                           border-radius: 10px;
                           }
                           .slider-input {
                           width: 100%;max-width: 700px;}
                          .slider-input .irs-bar, .slider-input .irs-bar-edge {
                          background-color: #007BFF; /* Blue slider bar */
                          border-color: #007BFF;}
                          .slider-input .irs-slider {
                          background-color: #ffffff;
                          border: 16px solid #007BFF;
                          }
                          .slider-input .irs-single {
                          background: #007BFF;
                          color: #ffffff;
                          font-weight: bold;
                          padding: 8px 12px;
                          border-radius: 4px;
                         /* Style for slider label */
                         .slider-label {
                         font-weight: bold;
                         font-size: 22px;}}"))
                         ),

                         fluidRow(
                           leafletOutput("map", height = 450, width = 700),
                           br(),

                           sliderInput("date_slider",
                                       label = tags$span("Select a Date:", class = "slider-label"),

                                       min = as.Date("2023-01-01"),
                                       max = as.Date("2023-12-31"),
                                       value = as.Date("2023-06-10"),
                                       step = 7,
                                       timeFormat = "%Y-%m-%d",
                                       animate = F,
                                       width = "100%"  # Full-width slider, max-width applied via CSS
                           ),
                           br(),
                           div(style = "text-align: center; font-weight: bold; color: #007BFF;",
                               ""),
                           br(),
                         ),

                     #    div(style = "display: flex; align-items: flex-start;",
                         #    materialSwitch(inputId = "checkbox_spectra", label = "Get Spectral Information", status = "danger"),
                          #   p('Activate it to see spectral profiles on the Spectral tab', style = "margin-left: 2px;")
                      #   ),

                ),
                tabPanel("Bandset",
                         plotOutput("plots_band"), # For trait distribution plot
                ),
                tabPanel("Spectral profile",

                        DTOutput("pixel_table"),
                        plotOutput("spectral_plot"),
                        br(),
                        p(strong('The selected satellite imagery used in the selected aggregation method:')),
                        br(),

                        DTOutput("satellite_table"),
                        # Row layout for download buttons
                        div(
                          style = "display: flex; align-items: center; gap: 10px; padding-top: 10px;",
                          downloadButton("download_table", "Download Spectral Data"),
                          downloadButton("download_metadata", "Download Metadata")
                        ),


                ),
                tabPanel("Info STAC",
                         br(),
                         p(style = "text-align: justify;", HTML("This Shiny app allows you to explore and analyze Sentinel-2 satellite imagery.
                         It takes advantage of cloud-optimized GeoTIFFs (COGs) and the SpatioTemporal Asset Catalog (STAC) for efficient data access.
                         The app can access Sentinel-2 catalogs from two providers.
                         <ul>
                         <li>Amazon Web Services: <a href='https://earth-search.aws.element84.com/v0/collections/sentinel-s2-l2a' target='_blank'>https://earth-search.aws.element84.com/v0/collections/sentinel-s2-l2a</a></li>
                         <li>Microsoft Planetary Computer: <a href='https://planetarycomputer.microsoft.com/api/stac/v1' target='_blank'>https://planetarycomputer.microsoft.com/api/stac/v1</a></li>
                         </ul>
                         This dual-source approach ensures seamless searching and retrieval of the imagery you need.")),

                         ),
                tabPanel("References",
                         p(style = "text-align: justify;",HTML('If you use <b>ToolsRTM</b> or <b>SCOPEinR</b> package, please cite the following references:')),

                         p(style = "text-align: justify;",HTML('Camino et al., (2024). RT-Simulator: An Online Platform to Simulate Canopy Reflectance from Biochemical and Structural Plant Properties Using Radiative Transfer Models,
                                   <i>IGARSS 2024 - 2024 IEEE International Geoscience and Remote Sensing Symposium</i>, Athens, Greece, 2024, pp. 2811-2814,
                                          <a href="https://ieeexplore.ieee.org/document/10642442" target="_blank">doi: 10.1109/IGARSS53475.2024.10642442</a>.')),

                         p(style = "text-align: justify;",HTML('Arano et al., (2024). Enhancing Chlorophyll Content Estimation With Sentinel-2 Imagery: A Fusion of Deep Learning and Biophysical Models,
                                   <i>IGARSS 2024 - 2024 IEEE International Geoscience and Remote Sensing Symposium</i>, Athens, Greece, 2024, pp. 4486-4489,
                                   <a href="https://ieeexplore.ieee.org/document/10641613" target="_blank">doi: 10.1109/IGARSS53475.2024.10641613</a>.')),

                         p('Camino et al., (in prep). Integrating physiological plant traits with Sentinel-2
                                     imagery for monitoring gross primary production and detecting forest disturbances. .'),
                )# end tabPanel
              )

    )
  )
 )
)

# 4. Define server logic  -----


server <- function(input, output,session) {

  # Get the paths of files in the "www/scenario/" directory using list.files()
  scenari_path <- (list.files("www/scenario/", full.names = TRUE))

  # Create a mapping from scenario names to file paths
  scenario_map <- list("Bark beetle outbreak" = scenari_path[1],
    "San Rossore" = scenari_path[2])


  print( scenari_path)
  # Reactive expression for selected scenario
  scenario <- reactive({
    req(input$scenario)  # Ensure input is available
    selected_path <- scenario_map[[input$scenario]]  # Get the path based on the selected scenario
    sf::read_sf(selected_path)  # Read the selected scenario
    # Calculate the bounding box using the buffer size around the scenario centroid

  })


  scenario_coord <- reactive({
    req(scenario())
    # Get the CRS of the scenario
    scenario_crs <- sf::st_crs(scenario())
    # Merge geometries if there are multiple
    merged_geom <- sf::st_union(scenario())

    coord_start <- merged_geom |>
      sf::st_centroid() |>
      sf::st_transform(4326) |>  # Ensuring it uses CRS 4326
      sf::st_coordinates()

    # Print coordinates and CRS to check
    print(paste("Coordinates:", paste(coord_start, collapse = ", ")))
    print(paste("CRS:", scenario_crs$proj4string))
    return(coord_start)
  })

  ## Leaflet map viewer -----------------------------------------------------
  output$map <- renderLeaflet({
    req(scenario())  # Ensure scenario is not NULL
    coord_start <- scenario_coord()  # Get the coordinates
    ## get the map

  # Create the Leaflet map
    m <- leaflet(data = scenario()) |>
      setView(lng = coord_start[1], lat = coord_start[2], zoom = 16) |>
      addTiles(urlTemplate = "https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}",
               group = "Orthophoto") |>  # Add Esri layer first
      addTiles(group = "OSM") |>          # Add OSM layer
      addPolylines(color = "#7df9ff", opacity = 1, weight = 3, layerId = 'scenario') |>
      addLayersControl(baseGroups = c("Orthophoto", "OSM"),  # Order here determines default visibility
                       options = layersControlOptions(collapsed = FALSE)) |>
      addScaleBar(position = "bottomright") |>
      addMeasure(primaryLengthUnit = "meters", primaryAreaUnit = "sqmeters")

    # Additional map customization based on selected spatial platform
    if (input$server_map == "Amazon Web Service") {

      m <- m |>
        addPopups(coord_start[1], coord_start[2],
                  "A spatial map using Amazon cloud services is display.",
                  options = popupOptions(closeButton = TRUE))


    } else if (input$server_map == "Microsoft Planetary Computer") {

      # Here, you can include any additional layers or customization specific to Microsoft Planetary Computer
      m <- m |>
        addPopups(coord_start[1], coord_start[2],
                  "Currently displaying data from Microsoft Planetary Computer.",
                  options = popupOptions(closeButton = TRUE))

    }
  })  |> bindEvent(input$scenario, input$server_map)


  # Reactive expression to get the collection ID based on user selection
  selected_collection <- reactive({
    req(input$satellite_collection)  # Ensure input is available

    # Get the date range and satellite collection from user inputs
   # date_range <- input$date_range_collection
    start_date <- input$date_slider -3
    end_date <- start_date + 7  # Add 10 days to the start date
    date_range <-c(start_date, end_date)
    satellite_collection <- input$satellite_collection
    cloud_cover <- input$cloud_collection
    # Map collection names to their respective Stack Catalog IDs


    if(input$server_map == 'Amazon Web Service'){
      collection_map <- list(
        "Sentinel-2 Collection" = 'sentinel-s2-l2a')
    } else if (input$server_map == 'Microsoft Planetary Computer'){
      collection_map <- list(
        "Sentinel-2 Collection" = 'sentinel-2-l2a')
    }

    # Get the corresponding Stack Catalog ID
    collection_id <- collection_map[[input$satellite_collection]]
    # Print the selected collection ID and date range to the console
    print(paste("Selected Collection ID:", collection_id))
    print(paste("Selected Date Range:", date_range))

    return(collection_id)

  })


  # Initialize reactive values to store raster objects
  avg_raster_val <- reactiveVal(NULL)
  ndvi_avg_val <- reactiveVal(NULL)



  observeEvent(input$get_cube, {
    show_modal_spinner()
    #req(input$date_range_collection)
    req(input$date_slider)
   # req(input$date_)
    req(scenario())  # Ensure scenario is not NULL
    req(scenario_coord())  # Get the coordinates

    print('scenario()')
    print(class(scenario()))

    scenario() |>
      sf::st_centroid() |>
      sf::st_transform(4326) |>  # Ensuring it uses CRS 4326
      sf::st_coordinates()
    print('coord_start')
    print(class(scenario()))


    # Get the selected collection ID from the reactive expression
    collection <- selected_collection()

    # Get the date range from user input

    start_date <- input$date_slider -3
    end_date <- start_date + 7  # Add 10 days to the start date
    date_range <-c(start_date, end_date)
    print(date_range)
    #date_range <- input$date_range_collection

    cloud_threshold <- input$cloud_collection
    print(paste0('cloud is:',cloud_threshold))

    showNotification("Getting the collection from Microsoft Planetary Computer ...", type = "message")


    if(input$server_map == 'Amazon Web Service'){
       cloud_ <- 'amazon'
    } else if (input$server_map == 'Microsoft Planetary Computer'){
      cloud_ <- 'microsoft'
    }

    print(collection)
    satellite_collection <- get.satellite_collection(scenario=scenario(), collection= collection,
                                                     cloud_server =cloud_,
                                                     date_range = date_range,
                                                     n.limit=3,
                                                     cloud_threshold= cloud_threshold,
                                                     buffer_size = input$bbox_cube)
    if(missing(satellite_collection)){
      showNotification("Error: Failed to retrieve the collection. Please check your inputs and try again.", type = "error")
      remove_modal_spinner()
    }
      # Check if satellite_collection is empty
    # Check if satellite_collection is NULL
    if (is.null(satellite_collection[[1]])) {
      showNotification("Error: Failed to retrieve the collection. Please check your inputs and try again.", type = "error")
      remove_modal_spinner()
      return()
    } else {
      # Print satellite collection for debugging
      print(head(satellite_collection[[2]]))
      showNotification("Collection retrieved successfully.", type = "message")

    }


    # Calculate the centroid of the scenario geometry
    centroid <- sf::st_centroid(scenario())

    # Apply a buffer around the centroid based on the specified buffer size
    buffer_area <- sf::st_buffer(centroid, input$bbox_cube)

    if (input$satellite_collection == 'Sentinel-2 Collection'){

      # Show notification that the process is starting
      showNotification("Fetching satellite cube for plotting. This may take some time...", type = "warning")
      # Show a progress bar while processing
      withProgress(message = 'Extracting satellite collection ...', {
        # Update progress
        total_steps <- 2
        incProgress(1/total_steps, detail = "\n")

        # Update choices based on selected method
        if (input$aggregation_method == "Minimum") {
          aggregation_ <- 'min'
          s_collection <- satellite_collection[[1]]

        } else if (input$aggregation_method == "Maximum") {
          aggregation_ <- 'max'
          s_collection <- satellite_collection[[1]]

        } else if (input$aggregation_method == "Mean") {
          aggregation_ <- 'mean'
          s_collection <- satellite_collection[[1]]

        } else if (input$aggregation_method == "Median") {
          aggregation_ <- 'median'
          s_collection <- satellite_collection[[1]]

        } else if (input$aggregation_method == "First") {
          aggregation_ <- 'first'
          s_collection <- satellite_collection[[1]]
        }

        print(aggregation_)
        # Update choices based on selected method
        if (input$resampling_method == "Near") {
          resampling_ <- 'near'
        } else if (input$resampling_method == "Bicubic") {
          resampling_ <- 'bicubic'
        } else if (input$resampling_method == "Bilinear") {
          resampling_ <- 'bilinear'
        }
        print(resampling_)


        avg_raster_cube <- get.sentinel2_cube(s_collection, shape = buffer_area,
                                              date_range = date_range,
                                              aggregation_method = aggregation_,
                                              resampling_method = resampling_,
                                              get.dataset = FALSE)
        print('raster done .....')
        # Final step
        incProgress(1/total_steps, detail = "Finalizing ...")
      })
      showNotification("Collection cube processed successfully.", type = "message")
      shinybusy::remove_modal_spinner()
      # Print summary and plot
      print(summary(avg_raster_cube))
      #plot(avg_raster_cube)

      # NB (raster->terra migration): get.sentinel2_cube(get.dataset = FALSE)
      # already returns a terra SpatRaster directly (see ToolsRTM::get.sentinel2_cube's
      # own source) -- the raster package's brick() was converting an already-terra object
      # into a legacy RasterBrick for no functional reason.
      avg_raster <- avg_raster_cube
      avg_raster_val <- avg_raster_val(avg_raster)
      # Calculate NDVI for average and median raster cubes
      showNotification("Calculating NDVI for raster layers ...", type = "message")
      ndvi_avg <- (avg_raster[[7]] - avg_raster[[3]]) / (avg_raster[[7]] + avg_raster[[3]])
      names(ndvi_avg) <- "NDVI"
      ndvi_avg_val <- ndvi_avg_val(ndvi_avg)


      # Notify the user that processing is complete
      showNotification("NDVI and raster layers processed successfully.", type = "message")
      # Assign names to NDVI rasters
    }

    # Render plot for average raster values using ggplot2
    output$plots_band <- renderPlot({
      req(avg_raster_cube)

      # Convert SpatRaster to data frame
      df <- as.data.frame(avg_raster_cube, xy = TRUE, na.rm = TRUE)

      if (input$satellite_collection == 'Sentinel-2 Collection' ){
        # Exclude the 'SCL' column, pivot longer, and multiply values by 1/10000
        df_long <- df %>%
          dplyr::select(-SCL) %>%  # Exclude the 'SCL' column
          pivot_longer(cols = starts_with("B"), names_to = "band", values_to = "value") %>%
          mutate(value = value * 1/10000)  # Multiply values by 1/10000

      }


      # Find the min and max values to set appropriate color scale limits
      min_value <- min(df_long$value, na.rm = TRUE)
      max_value <- max(df_long$value, na.rm = TRUE)

      # Plot with adjusted color scale limits
      ggplot(df_long, aes(x = x, y = y, fill = value)) +
        geom_tile() +
        coord_equal() +
        viridis:::scale_fill_viridis(limits = c(min_value, max_value), na.value = "gray90") +   # Set limits and color for NA values
        facet_wrap(~ band, ncol = 5) +  # Facet by band
        labs(title = "", fill = "Reflectance") +
        theme_bw() +
        theme(
          plot.title = element_text(hjust = 0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),  # Remove x-axis text
          axis.text.y = element_blank(),  # Remove y-axis text
          axis.ticks = element_blank() ,    # Remove axis ticks
          strip.text = element_text(size = 16),  # Increase facet label size to 16
          legend.text = element_text(size = 16)   # Increase legend text size to 16
        )
    })

    # Set the raster_ready output to TRUE to show the download button
    output$raster_ready <- reactive({ TRUE })
    outputOptions(output, "raster_ready", suspendWhenHidden = FALSE)  # Ensure reactivity

    # Now render the Leaflet map with the computed avg_raster_cube
    output$map <- renderLeaflet({
      req(scenario())  # Ensure scenario is not NULL
      coord_start <- scenario_coord()  # Get the coordinates

      # Initialize the map with scenario and OSM layer
      m <- leaflet(data = scenario()) |>
        clearControls() |>
        clearGroup(c("Orthophoto", "OSM", "RGB", "False Color", "NDVI")) |> # Clear previous groups
        setView(lng = coord_start[1], lat = coord_start[2], zoom = 16) |>
        clearControls() |>  # Clears previous legends and controls
        addTiles(urlTemplate = "https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}",
                 group = "Orthophoto") |>  # Add Esri layer first
        addTiles(group = "OSM") |>          # Add OSM layer
        addPolylines(color = "#7df9ff", opacity = 1, weight = 3, layerId = 'scenario') |>
        addLayersControl(baseGroups = c("Orthophoto", "OSM"),  # Order here determines default visibility
                         options = layersControlOptions(collapsed = FALSE)) |>
        addScaleBar(position = "bottomright") |>
        addMeasure(primaryLengthUnit = "meters", primaryAreaUnit = "sqmeters")


      # Add the average raster RGB layer (after it's computed)
      if (!is.null(avg_raster_cube)) {
        # Define a viridis color palette for NDVI
        ndvi_pal <- colorNumeric(palette = "viridis", domain = c(-0.3, 0.9),  reverse = TRUE)

        m <- m |>
          leafem::addRasterRGB(avg_raster, r = 1, g = 2, b = 3,
                               quantile = c(0.02, 0.98),
                               group = "RGB") |>
          leafem::addRasterRGB(avg_raster, r = 7, g = 4, b = 3,
                               quantile = c(0.1, 0.9),
                               group = "False Color") |>
          # Add NDVI layers
          #leaflet::addRasterImage(ndvi_avg, group = "NDVI", colors = terrain.colors(100))
          leaflet::addRasterImage(ndvi_avg, group = "NDVI", colors = ndvi_pal) %>%
          leaflet::addLegend(pal = ndvi_pal, values = c(-0.3,0.9), title = "NDVI")


      }

      # Add layer controls for the raster and base layers
      m <- m |>

        addLayersControl(
          baseGroups = c("Orthophoto", "OSM"),
          overlayGroups = c("RGB", "False Color", "NDVI"),  # Add other layers as needed
          options = layersControlOptions(collapsed = FALSE)
        ) |>
        addScaleBar(position = "bottomright") |>
        addMeasure(primaryLengthUnit = "meters", primaryAreaUnit = "sqmeters") |>
        leafem::addMouseCoordinates()

      # Add click event to extract pixel values
      m <- m |>
        htmlwidgets::onRender("
        function(el, x) {
          this.on('click', function(e) {
            var lat = e.latlng.lat;
            var lng = e.latlng.lng;

            // Send coordinates back to Shiny
            Shiny.setInputValue('click_coordinates', {lat: lat, lng: lng});

            // Display coordinates in a popup
            L.popup()
              .setLatLng(e.latlng)
              .setContent('Clicked coordinates: <br>Latitude: ' + lat + '<br>Longitude: ' + lng)
              .openOn(this);
          });
        }"
        )

      # Return the map
      m
    })

    #shinybusy::remove_modal_spinner()


    # In your server function
    pixel_data <- reactiveVal(data.frame(Latitude = numeric(), Longitude = numeric(), PixelInfo = character()))

    leaflet_map <- reactiveVal(NULL)
    # Initialize a counter for marker IDs
    marker_counter <- reactiveVal(0)

    #Add this observeEvent somewhere in your server code
    observeEvent(input$click_coordinates, {
      lat <- input$click_coordinates$lat
      lng <- input$click_coordinates$lng
      # Print coordinates to the R console
      cat("Clicked coordinates: Latitude:", lat, "Longitude:", lng, "\n")
      # Increment the counter and generate the marker ID
      marker_counter(marker_counter() + 1)
      marker_id <- paste0("", marker_counter())

      # Update the map with a marker with the ID
      leafletProxy("map") %>%  # Use leafletProxy directly
        addMarkers(lng = lng, lat = lat,
                   popup = paste("ID:", marker_id, "<br>",  # Include ID in the popup
                                 "Latitude:", lat, "<br>", "Longitude:", lng),
                   layerId = marker_id,      label = lapply(marker_id, htmltools::HTML))

      # Show notification to user
      showNotification("Spectral information has been saved in the Spectral tab.",
                       type = "message", duration = 3)
      # Create a data frame for the extraction coordinates
      extraction_coords <- data.frame(Latitude = lat, Longitude = lng)

      # Get pre-calculated raster layers
      avg_raster <- avg_raster_val()  # Average raster data
      ndvi_raster <- ndvi_avg_val()    # NDVI raster data

      # Create a data frame for the extraction coordinates
      # Add this line to convert the clicked coordinates to the same CRS as the raster
      extraction_coords_sf <- st_as_sf(extraction_coords, coords = c("Longitude", "Latitude"), crs = 4326)  # Assuming WGS84 (EPSG:4326)
      extraction_coords_transformed <- st_transform(extraction_coords_sf, st_crs(avg_raster))

      # Extract pixel information using terra::extract
      # Now extract using transformed coordinates
      raster_values <- terra::extract(avg_raster, st_coordinates(extraction_coords_transformed))
      ndvi_values <- terra::extract(ndvi_raster, st_coordinates(extraction_coords_transformed))

      print(raster_values)

      # Print extracted values for debugging
      cat("Extracted raster values:", raster_values, "\n")
      cat("Extracted NDVI values:", ndvi_values, "\n")

      # Check if extraction returned any values
      if (is.null(raster_values) || is.null(ndvi_values)) {
        cat("No data found for the clicked coordinates.\n")
        return()
      }

      # Combine extracted data into a single data frame
      factor_rfl = 1/10000
      # Create a data frame for pixel information
      pixel_info <- data.frame(
        ID = marker_id,
        Latitude = lat,
        Longitude = lng,
        B02 = raster_values[1, "B02"] * factor_rfl,
        B03 = raster_values[1, "B03"] * factor_rfl,
        B04 = raster_values[1, "B04"] * factor_rfl,
        B05 = raster_values[1, "B05"] * factor_rfl,
        B06 = raster_values[1, "B06"] * factor_rfl,
        B07 = raster_values[1, "B07"] * factor_rfl,
        B08 = raster_values[1, "B08"] * factor_rfl,
        B8A = raster_values[1, "B8A"] * factor_rfl,
        B11 = raster_values[1, "B11"] * factor_rfl,
        B12 = raster_values[1, "B12"] * factor_rfl,
        SCL = raster_values[1, "SCL"],
        NDVI = ndvi_values  # Assuming ndvi_values is a single value
      )

      # Print pixel info for debugging
      cat("Pixel Info:\n")
      rownames(pixel_info) <-NULL
      # Save the new data to the reactive value
      current_data <- pixel_data()
      updated_data <- rbind(current_data, pixel_info)
      pixel_data(updated_data)
      print(pixel_data)
    })

    # get the spectral table
    output$pixel_table <- renderDT({
     # req(input$checkbox_spectra)  # Only render if the switch is active
      # Get the pixel data
      pixel_data()
    })

    # get the metadata table
    output$satellite_table <- renderDT({
     # req(input$checkbox_spectra)  # Only render if the switch is active
      # Get the pixel data
      satellite_collection[[2]]
    })

    # Download the data table
    output$download_table <- downloadHandler(
      filename = function() {
        "pixel_data.csv"
      },
      content = function(file) {
        write.csv(pixel_data(), file, row.names = FALSE)
      }
    )

    # Download the data table
    output$download_metadata <- downloadHandler(
      filename = function() {
        "metadata.csv"
      },
      content = function(file) {
        df.metadata <-  satellite_collection[[2]]
        write.csv(df.metadata, file, row.names = FALSE)
      }
    )

    # Reactive value to store pixel data
    pixel_data <- reactiveVal(data.frame())

    # Render the spectral plot
    output$spectral_plot <- renderPlot({
      #req(input$checkbox_spectra, pixel_data()) # Require switch to be active and data available

      # Prepare data for plotting
      plot_data <- pixel_data() %>%
        tidyr::pivot_longer(cols = c("B02", "B03", "B04", "B05", "B06", "B07", "B08", "B11", "B12"),
                            names_to = "Band", values_to = "Reflectance")

      # Create the ggplot
      ggplot(plot_data, aes(x = Band, y = Reflectance, group = ID, color = factor(ID))) +
        geom_line() +   geom_point() +
        labs(title = "", x = "", y = "Reflectance", color = "ID") +  # Add legend title
        theme_bw() + theme(
          legend.position = 'right',
          strip.text = element_text(face = "bold", size = 14),  # Bold facet labels
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.title = element_text(face = "bold", size = 14),
          axis.text.y = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x = element_text(hjust = 0.5, size = 14, face = "bold"),
          panel.grid.major = element_blank(),  # Optional: Remove grid lines
          panel.grid.minor = element_blank()) +

        scale_color_discrete(name = "ID")
    })

      # Observe changes in scenarios
  observeEvent(input$scenario, {
    # Reset pixel_data
    pixel_data(data.frame())

    # Reset marker counter
    marker_counter(0)

    # Clear markers from the map
    leafletProxy("map") %>% clearMarkers()
  })
    # Enable download of the avg_raster_cube as a GeoTIFF
    output$download_tiff <- downloadHandler(
      filename = function() {
        paste("satellite_product", Sys.Date(), ".tiff", sep = "")
      },
      content = function(file) {
        req(avg_raster_val())  # Ensure avg_raster is not NULL
        avg_raster_val <- avg_raster_val()  # Get the pre-calculated avg_raster
        ndvi_avg_val <- ndvi_avg_val()  # Get the pre-calculated ndvi_avg

        # NB (raster->terra migration): avg_raster/ndvi_avg are terra
        # SpatRasters now (avg_raster_cube comes straight from
        # get.sentinel2_cube(); ndvi_avg is SpatRaster arithmetic on it) --
        # the old inherits(ndvi_avg, "RasterLayer") check was always FALSE
        # post-migration, silently taking the "NDVI not available" branch
        # every time. Checking for "SpatRaster" instead restores the
        # intended behavior.
        if (!is.null(ndvi_avg) && inherits(ndvi_avg, "SpatRaster")) {
          # If NDVI is available, combine avg_raster and ndvi_avg into one
          # multi-layer SpatRaster (terra's equivalent of the raster package's stack()).
          combined_raster <- c(avg_raster, ndvi_avg)

          # Optionally, you can set names for the combined raster layers
          names(combined_raster) <- c(names(avg_raster), "NDVI")

          # Write the combined raster to a GeoTIFF file
          terra::writeRaster(combined_raster, file, filetype = "GTiff", overwrite = TRUE)

        } else {
          # If NDVI is not available, write only the avg_raster
          terra::writeRaster(avg_raster, file, filetype = "GTiff", overwrite = TRUE)
        }
      }
    )

    observe({
      if (input$scenario == 'Bark beetle outbreak' | input$scenario == 'San Rossore') {
        showNotification("Scenario loaded successfully.", type = "message")
      }


    })




  })



}
# 5. Create Shiny object-----

shinyApp(ui = ui, server = server)

