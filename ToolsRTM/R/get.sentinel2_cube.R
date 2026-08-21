#' Create and Process a Sentinel-2 Data Cube
#'
#' This function generates a data cube from a Sentinel-2 image collection by aggregating the pixel values over a specified spatial extent and date range.
#'
#' @param s_collection An object representing the Sentinel-2 image collection, typically created using a function like `get.satellite_collection()`.
#' @param shape An object defining the spatial extent for the data cube, such as a polygon or a raster extent.
#' @param date_range A character vector of length 2 specifying the start and end dates in "YYYY-MM-DD" format for the aggregation period.
#' @param aggregation_method A character string indicating the method of aggregation to apply to the pixel values. Default is "mean". Other options may include "min", "max", "median", or "first"
#' @param resampling_method A character string indicating the method of resampling to apply to the pixel values. Default is "bicubic". Other options may include "near", or "bilinear".
#' @param get.dataset A logical value indicating whether to return the dataset of processed data (default is TRUE).
#'
#' @return A processed data cube containing aggregated pixel values over the specified shape and date range. If `get.dataset` is TRUE, it returns a dataset object; otherwise, it returns the cube object.
#' @export
#'
#' @examples
#' \dontrun{
#' # Example usage of the get.sentinel2_cube function
#'
#' # Assuming s_collection has been obtained from a previous call
#' s_collection <- get.satellite_collection("your_scenario_here", "sentinel-2-l2a",
#'                                           c("2023-01-01", "2023-01-31"),
#'                                           cloud_threshold = 20)
#'
#' # Define the spatial shape (e.g., a polygon)
#' poly <- list(rbind(c(4.8, 52.2), c(4.8, 52.4), c(5.0, 52.4), c(5.0, 52.2), c(4.8, 52.2)))
#' shape <- sf::st_as_sf(data.frame(id = 1,
#'                                   geometry = sf::st_sfc(sf::st_polygon(poly))),
#'                        crs = 4326)
#'
#' # Define date range for aggregation
#' date_range <- c("2023-01-01", "2023-01-31")
#'
#' # Retrieve the Sentinel-2 data cube
#' sentinel2_cube <- get.sentinel2_cube(s_collection, shape, date_range,
#'                                       aggregation_method = "mean", get.dataset = TRUE)
#'
#' # Print the resulting data cube details
#' print(sentinel2_cube)
#' }
get.sentinel2_cube <- function(s_collection, shape, date_range, 
                               aggregation_method = "mean", 
                               resampling_method='bicubic', get.dataset=T) {

  
  if (missing(get.dataset) || is.null(get.dataset)) {
    get.dataset = F  # Set default to F if empty or missing
    print('A raster is processing ---')
  } else {
    print('The reflectance are processing ---')
  }

  crs_cube <- "EPSG:3035"  # Use a projected CRS for raster operations


  # Ensure the shape is in sf format
  if (!inherits(shape, "sf")) {
    shape <- sf::st_as_sf(shape)
  }
  # Convert bounding box to sf object and transform to the target CRS
  shape_cube <- shape |>
    sf::st_as_sfc() |>
    sf::st_transform(crs_cube) |>
    sf::st_bbox(crs = crs_cube)


  # Set default values if missing or NULL
  if (missing(aggregation_method) || is.null(aggregation_method)) {
    aggregation_method <- "mean"
  }
  
  if (missing(resampling_method) || is.null(resampling_method)) {
    resampling_method <- "bicubic"
  }
  
  # Define the cube view (spatial extent and resolution)
  view <- gdalcubes::cube_view(
    srs = crs_cube,
    extent = list(
      t0 = as.character(date_range[[1]]),
      t1 = as.character(date_range[[2]]),
      left = shape_cube["xmin"],
      right = shape_cube["xmax"],
      top = shape_cube["ymax"],
      bottom = shape_cube["ymin"]
    ),
    dx = 20, dy = 20, dt = "P1D",  # 20m resolution, time unit = 1 day
    aggregation = aggregation_method, #"aggregation method
    resampling = resampling_method    # Resampling method
  )


  # Mask to remove clouds and shadows using the SCL band
  mask <- gdalcubes::image_mask("SCL", values = c(0,1,3,8,9,10,11,12))  # Clouds,snow,water bad pixels and shadows

  .gc_par <- max(1, parallel::detectCores()-2)
  if (nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_"))) .gc_par <- min(.gc_par, 2L)
  gdalcubes::gdalcubes_options(parallel = .gc_par)
  # Create a raster cube for the given collection and view
  cube <- gdalcubes::raster_cube(s_collection, view,mask)
  
  # Check available bands in the cube
  available_bands <- names(cube)  # Adjust based on how you access bands

  # Define the preferred order of bands
  preferred_order <- c("B02", "B03", "B04", "B05", "B06", "B07", "B08", "B8A","B11","B12",'SCL')
  # Reorder the available bands based on preferred order
  selected_bands <- preferred_order[preferred_order %in% available_bands]
  print(selected_bands)

  # Calculate average or median across time for selected bands
  # Calculate the result cube based on the selected aggregation method
  if (aggregation_method == "mean") {
    result_cube <- cube |>
      gdalcubes::select_bands(selected_bands) |>
      gdalcubes::reduce_time(paste0("mean(", selected_bands, ")"))
    
  } else if (aggregation_method == "median") {
    result_cube <- cube |>
      gdalcubes::select_bands(selected_bands) |>
      gdalcubes::reduce_time(paste0("median(", selected_bands, ")"))
    
  } else if (aggregation_method == "min") {
    result_cube <- cube |>
      gdalcubes::select_bands(selected_bands) |>
      gdalcubes::reduce_time(paste0("min(", selected_bands, ")"))
    
  } else if (aggregation_method == "max") {
    result_cube <- cube |>
      gdalcubes::select_bands(selected_bands) |>
      gdalcubes::reduce_time(paste0("max(", selected_bands, ")"))
    
  } else if (aggregation_method == "first") {
    result_cube <- cube |>
      gdalcubes::select_bands(selected_bands) |>
      #gdalcubes::reduce_time(paste0("min(", selected_bands, ")"))
      gdalcubes::slice_time(as.character(date_range[[1]]))
    
  } else {
    # Default case if aggregation_method is invalid
    warning("Invalid aggregation method specified. Using 'mean' as the default.")
    result_cube <- cube |>
      gdalcubes::select_bands(selected_bands) |>
      gdalcubes::reduce_time(paste0("mean(", selected_bands, ")"))
  }
  
 

  # Convert to stars or terra format for visualization
  # NB (fix): bare rast() only resolves if the CALLER happens to have
  # library(terra) attached (not just installed) -- terra is Imported by
  # this package but rast() was never called with its namespace prefix, so
  # any caller with only `requireNamespace("terra")`/no attach hit
  # "could not find function \"rast\"" even though the package itself
  # declares and uses terra internally elsewhere (e.g. get.satellite_collection()).
  result_raster <- gdalcubes::st_as_stars.cube(result_cube) |> terra::rast()

  # Generate names based on bands and time range
  band_names <- paste0(selected_bands, "_",
                       format(as.Date(date_range[[1]]), "%Y%m%d"), "_",
                       format(as.Date(date_range[[2]]), "%Y%m%d"))

  # Assign names to the raster layers
  names(result_raster) <- selected_bands

  if (get.dataset == T){
    data.collection <- cube |> as.data.frame(complete_only = TRUE)
    data.list <-list(data.collection=data.collection, raster=result_raster)

    return(list(result_raster,data.list))

  } else {
    return(result_raster)
  }


}
