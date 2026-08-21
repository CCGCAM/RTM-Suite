#' Retrieve Satellite Collections
#'
#' This function retrieves satellite image collections from the Planetary Computer STAC API based on specified parameters, including bounding box, date range, and cloud cover threshold.
#'
#' @param scenario A character string representing the scenario, which is used to obtain the bounding box.
#' @param collection A character string specifying the name of the satellite data collection (e.g., "sentinel-2-l2a").
#' @param collection A character string specifying the name of the satellite data collection (e.g., "sentinel-2-l2a").
#' @param cloud_server The cloud server to use. Options are "microsoft" or "amazon". Defaults to "microsoft".
#' @param n.limit a integer with the maximum  number of images during the extractions.
#' @param cloud_threshold A numeric value representing the maximum allowable cloud cover percentage.
#' @param buffer_size A numeric value indicating the buffer size around the bounding box in meters (default is 300).
#'
#' @return A filtered image collection based on the specified criteria.
#' @export
#'
#' @examples
#' \dontrun{
#' # Example usage of the get.satellite_collection function
#' scenario <- "your_scenario_here"  # Replace with your actual scenario
#' # names ofr every single cloud server
#' coleccion.microsoft <- c('sentinel-2-l2a','landsat-c2-l2','modis-17A2HGF-061','modis-09A1-061','modis-09Q1-061','modis-11A2-061','modis-15A2H-061','modis-15A3H-061')
#' coleccion.aws <- c('sentinel-s2-l2a','sentinel-s2-l2a-cogs','landsat-8-l1-c1')
#' collection <- "sentinel-2-l2a"  # Sentinel-2 Level-2A collection
#' date_range <- c("2023-06-01", "2023-09-30")  # Date range for January 2023
#' cloud_threshold <- 20  # Maximum cloud cover of 20%
#'
#' # Retrieve the satellite collection
#' satellite_data <- get.satellite_collection(scenario, collection, date_range, cloud_threshold)
#'
#' # Print the retrieved collection details
#' print(satellite_data)
#' }
#'
get.satellite_collection <- function(scenario, collection, cloud_server='microsoft', n.limit=NULL,
                                     date_range, cloud_threshold, buffer_size = NULL) {
  
  if (missing(n.limit)) {
    n.limit = 500
  }
  gdalcubes::gdalcubes_options(parallel = parallel::detectCores()-2)
  # Ensure the scenario is in sf format
  if (!inherits(scenario, "sf")) {
    scenario <- sf::st_as_sf(scenario)
  }
  # Get bounding box based on the scenario centroid with a buffer of buffer_size meters
  bbox <- get_bounding_box(scenario, buffer_size)


  if (missing(cloud_server) || cloud_server == "microsoft") {

    stac.api <- 'Microsoft'
    planetary_computer <- rstac::stac("https://planetarycomputer.microsoft.com/api/stac/v1")
    message("STAC API from Microsoft is selected.")

  } else if (cloud_server == "amazon") {

    stac.api <- 'Amazon'
    planetary_computer <- rstac::stac("https://earth-search.aws.element84.com/v0")
    message("STAC API from Amazon is selected")

  } else {
    stop("Invalid cloud_server value. Please choose 'microsoft' or 'amazon'.")
  }


  if (cloud_server == "amazon"){

    # Search for items using the planetary computer STAC API
    items <- planetary_computer |>
      rstac::stac_search(collections = collection,
                  bbox = bbox[1:4],  # Use the reduced bounding box
                  datetime = paste0(date_range[[1]], "/", date_range[[2]]),
                  limit = n.limit) |>
      rstac::post_request() |> rstac::items_fetch(progress = FALSE)
      print(length(items$features))

  } else {
    # Search for items using the planetary computer STAC API
    items <- planetary_computer |>
      rstac::stac_search(collections = collection,
                  bbox = bbox[1:4],  # Use the reduced bounding box
                  datetime = paste0(date_range[[1]], "/", date_range[[2]]),
                  limit = n.limit) |>
      rstac::post_request()
    print(length(items$features))
  }

  # NB (fix): stac_search()/post_request()/items_fetch()/items_sign()/
  # sign_planetary_computer()/stac_image_collection() were all called bare
  # (no rstac::/gdalcubes:: prefix) throughout this function -- they only
  # ever resolved if the CALLER happened to have both library(rstac) and
  # library(gdalcubes) attached first (neither package is even in this
  # package's Imports, only Suggests), so a fresh call to
  # get.satellite_collection() with just library(ToolsRTM) attached failed
  # with e.g. "could not find function \"post_request\"".
  # Sign the STAC items for authentication
  items <- rstac::items_sign(items,
    rstac::sign_planetary_computer())

  
  if (length(items$features) == 0) {
    return(list(clt = NULL, df.data = NULL))
    stop("No images found in the specified collection for the given bounding box and date range.")
  } else {
    cat("Images found:", length(items$features), "\n")
  }

  # Get the assets for the first item
  first_item_assets <- items$features[[1]]$assets
  asset_names <- names(first_item_assets)
  print("Asset names:")  # Debug: Check the available asset names
  print(asset_names)

  # Define asset names based on the collection
  assets <- switch(collection,
                   ## AWS
                   "sentinel-s2-l2a" = c('B02','B03','B04','B05','B06','B07','B08','B8A','B11','B12','SCL'), # Example for Sentinel-2 Level-2A
                   "sentinel-s2-l2a-cogs" = c('B02','B03','B04','B05','B06','B07','B08','B8A','B11','B12','SCL'), # Example for Sentinel-2 Level-2A
                   ## Microsoft
                   "sentinel-2-l2a" = c('B02','B03','B04','B05','B06','B07','B08','B8A','B11','B12','SCL'), # Example for Sentinel-2 Level-2A
                   "landsat-c2-l2" = c('blue','green','red','nir08','swir16','swir22'), # for landsat level-2
                   "modis-17A2HGF-061" = c("Gpp_500m",'PsnNet_500m'), # Example for MODIS Gross Primary Productivity
                   "modis-11A2-061" = c("LST_Day_1km"), # or MODIS Land Surface Temperature
                   "modis-09A1-061" = asset_names[c(3:9)], #for MODIS Surface Reflectance 8-Day (500m)
                   "modis-09Q1-061" = c( "sur_refl_b01","sur_refl_b02","sur_refl_qc_250m","sur_refl_state_250m"), # Example for MODIS Surface Reflectance 8-Day (250m)
                   "modis-15A2H-061" = c("Lai_500m","Fpar_500m","LaiStdDev_500m",'FparStdDev_500m','FparLai_QC'), # Example for MODIS Leaf Area Index/FPAR 8-Day
                   "modis-15A3H-061" = c("Lai_500m","Fpar_500m","LaiStdDev_500m",'FparStdDev_500m','FparLai_QC')  # Example for MODIS Leaf Area Index/FPAR 4-Day
  )

  # Determine if cloud filtering is applicable
  if (collection == "modis-09A1-061" || collection == "modis-17A2HGF-061" || collection == "modis-11A2-061"
      || collection == 'modis-09Q1-061'|| collection == 'modis-09A1-061'
      || collection == 'modis-15A2H-061'|| collection == 'modis-15A3H-061') {
    # Create image collection without cloud cover filtering
  
    
    clt <- tryCatch(
      items$features |>
        gdalcubes::stac_image_collection(asset_names = assets),
      error = function(e) {
        NULL  # Return NA if an error occurs
      }
    )
    
    
  } else if (collection == 'landsat-c2-l2' || 
             collection == "sentinel-s2-l2a" || 
             collection == "sentinel-s2-l2a-cogs" || 
             collection == "sentinel-2-l2a") {
    # Create image collection with cloud cover filtering for other collections
    # Create image collection with cloud cover filtering for collections that require it
    clt <- tryCatch(
      items$features |>
        gdalcubes::stac_image_collection(asset_names = assets,
                              property_filter = function(x) {
                                x[["eo:cloud_cover"]] < cloud_threshold
                              }),
      error = function(e) {
        NULL  # Return NA if an error occurs
      }
    )
    
    
  } else {
    clt <- tryCatch(
      items$features |>
        gdalcubes::stac_image_collection(asset_names = assets),
      error = function(e) {
        NA  # Return NA if an error occurs
      }
    )
    
  }
  
  if (items$features[[1]]$collection == 'sentinel-2-l2a'){
    # If you want to extract the 'datetime' for further use, you can do that separately:
    id <- sapply(items$features, function(x) x$id)
    datetimes <- sapply(items$features, function(x) x$properties$datetime)
    platform <- sapply(items$features, function(x) x$properties$platform)
    zenith <- sapply(items$features, function(x) x$properties$`s2:mean_solar_zenith`)
    azimth <- sapply(items$features, function(x) x$properties$`s2:mean_solar_azimuth`)
    df.data <-data.frame(id=id, datetimes=datetimes, platform=platform,zenith=zenith,azimth=azimth)
    
  }  else if (items$features[[1]]$collection == 'sentinel-s2-l2a' || 
              items$features[[1]]$collection == 'sentinel-s2-l2a-cogs'){
    # If you want to extract the 'datetime' for further use, you can do that separately:
    id <- sapply(items$features, function(x) x$id)
    datetimes <- sapply(items$features, function(x) x$properties$datetime)
    platform <- sapply(items$features, function(x) x$properties$platform)
    gsd <- sapply(items$features, function(x) x$properties$gsd)
    off_nadir <- sapply(items$features, function(x) x$properties$`view:off_nadir`)
    df.data <-data.frame(id=id,datetimes=datetimes, platform=platform,off_nadir=off_nadir,gsd=gsd)
    
  } else if (items$features[[1]]$collection == 'landsat-c2-l2'){
    # If you want to extract the 'datetime' for further use, you can do that separately:
    id <- sapply(items$features, function(x) x$id)
    datetimes <- sapply(items$features, function(x) x$properties$datetime)
    platform <- sapply(items$features, function(x) x$properties$platform)
    zenith <- sapply(items$features, function(x) x$properties$`view:sun_azimuth`)
    azimth <- sapply(items$features, function(x) x$properties$`view:sun_azimuth`)
    df.data <-data.frame(id=id,datetimes=datetimes, platform=platform,zenith=zenith,azimth=azimth)
    
  } else if (items$features[[1]]$collection == 'modis-17A2HGF-061' || 
             items$features[[1]]$collection == 'modis-09A1-061'|| 
             items$features[[1]]$collection == "modis-11A2-061"|| 
             items$features[[1]]$collection == 'modis-09Q1-061'|| 
             items$features[[1]]$collection == 'modis-15A2H-061'|| 
             items$features[[1]]$collection == 'modis-15A3H-061'){
    # If you want to extract the 'datetime' for further use, you can do that separately:
    id <- sapply(items$features, function(x) x$id)
    datetimes <- sapply(items$features, function(x) x$properties$updated)
    platform <- sapply(items$features, function(x) x$properties$platform)
    df.data <-data.frame(id=id,datetimes=datetimes, platform=platform)
    
    
  }
  
  if (items$features[[1]]$collection == 'sentinel-2-l2a' || items$features[[1]]$collection == 'landsat-c2-l2' || 
      items$features[[1]]$collection == 'sentinel-s2-l2a' || items$features[[1]]$collection == 'sentinel-s2-l2a-cogs'){
    # Example: Select a specific date from the collection (e.g., the first image)
    selected_date <- as.character(df.data$datetimes[1])  # Change the index as needed
    selected_id <- as.character(df.data$id[1])  # Change the index as needed
    
    # Create the filtered collection for cloud cover and the specific date
    clt_first <- items$features |>
      gdalcubes::stac_image_collection(
        asset_names = assets,
        property_filter = function(x) {
          # Check if the cloud cover is below the threshold
          selected_id <- x[["id"]] == selected_id
          cloud_condition <- x[["eo:cloud_cover"]] < cloud_threshold
          # Check if the datetime matches the selected date
          date_condition <- x[["datetime"]] == selected_date
        }) 
    
  
  } else {
    # Example: Select a specific date from the collection (e.g., the first image)
    selected_date <- as.character(df.data$datetimes[1])  # Change the index as needed
    selected_id <- as.character(df.data$id[1])  # Change the index as needed
    
    # Create the filtered collection for cloud cover and the specific date
    clt_first <- items$features |>
      gdalcubes::stac_image_collection(
        asset_names = assets,
        property_filter = function(x) {
          # Check if the cloud cover is below the threshold
          selected_id <- x[["id"]] == selected_id
          # Check if the datetime matches the selected date
          date_condition <- x[["updated"]] == selected_date
        }) 
    
  }
  

  #
  return(list(clt,df.data,clt_first))  # Return the filtered image collection
}
