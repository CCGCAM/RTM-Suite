#' Define a Bounding Box Around a Scenario Centroid
#'
#' This function calculates a bounding box around the centroid of a given scenario geometry, applying a specified buffer size to determine the area of interest.
#'
#' @param scenario An object representing the scenario geometry (e.g., an `sf` object). The centroid of this geometry will be calculated to define the bounding box.
#' @param buffer_size A numeric value indicating the size of the buffer (in meters) to be applied around the centroid.
#'
#' @return A numeric vector of length 4 representing the bounding box in the format (xmin, ymin, xmax, ymax), projected to EPSG:4326 (WGS 84) coordinate system for compatibility with the STAC API.
#' @export
#'
#' @examples
#' \dontrun{
#' # Example usage of the get_bounding_box function
#' library(sf)
#'
#' # Create an example scenario geometry (a point in this case)
#' scenario <- st_point(c(4.9, 52.3)) |> st_sfc(crs = 4326)
#'
#' # Define buffer size in meters
#' buffer_size <- 500  # 500 meters
#'
#' # Retrieve the bounding box
#' bounding_box <- get_bounding_box(scenario, buffer_size)
#'
#' # Print the bounding box
#' print(bounding_box)
#' }
#'
get_bounding_box <- function(scenario, buffer_size) {

  # Ensure the scenario is in sf format
  if (!inherits(scenario, "sf")) {
    scenario <- sf::st_as_sf(scenario)
  }
  # Calculate the centroid of the scenario geometry
  centroid <- sf::st_centroid(scenario)

  # Apply a buffer around the centroid based on the specified buffer size
  buffer_area <- sf::st_buffer(centroid, buffer_size)

  # Get the bounding box of the buffered area
  bb <- sf::st_bbox(buffer_area)

  # Reproject the bounding box to EPSG:4326 (WGS 84) for STAC API compatibility
  bb_4326 <- sf::st_as_sfc(bb) |>
    sf::st_transform(4326) |>
    sf::st_bbox()

  return(bb_4326)  # Return the reprojected bounding box
}
