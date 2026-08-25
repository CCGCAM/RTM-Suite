# Create and Process a Sentinel-2 Data Cube

This function generates a data cube from a Sentinel-2 image collection
by aggregating the pixel values over a specified spatial extent and date
range.

## Usage

``` r
get.sentinel2_cube(
  s_collection,
  shape,
  date_range,
  aggregation_method = "mean",
  resampling_method = "bicubic",
  get.dataset = T
)
```

## Arguments

- s_collection:

  An object representing the Sentinel-2 image collection, typically
  created using a function like
  [`get.satellite_collection()`](get.satellite_collection.md).

- shape:

  An object defining the spatial extent for the data cube, such as a
  polygon or a raster extent.

- date_range:

  A character vector of length 2 specifying the start and end dates in
  "YYYY-MM-DD" format for the aggregation period.

- aggregation_method:

  A character string indicating the method of aggregation to apply to
  the pixel values. Default is "mean". Other options may include "min",
  "max", "median", or "first"

- resampling_method:

  A character string indicating the method of resampling to apply to the
  pixel values. Default is "bicubic". Other options may include "near",
  or "bilinear".

- get.dataset:

  A logical value indicating whether to return the dataset of processed
  data (default is TRUE).

## Value

A processed data cube containing aggregated pixel values over the
specified shape and date range. If `get.dataset` is TRUE, it returns a
dataset object; otherwise, it returns the cube object.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example usage of the get.sentinel2_cube function

# Assuming s_collection has been obtained from a previous call
s_collection <- get.satellite_collection("your_scenario_here", "sentinel-2-l2a",
                                          c("2023-01-01", "2023-01-31"),
                                          cloud_threshold = 20)

# Define the spatial shape (e.g., a polygon)
poly <- list(rbind(c(4.8, 52.2), c(4.8, 52.4), c(5.0, 52.4), c(5.0, 52.2), c(4.8, 52.2)))
shape <- sf::st_as_sf(data.frame(id = 1,
                                  geometry = sf::st_sfc(sf::st_polygon(poly))),
                       crs = 4326)

# Define date range for aggregation
date_range <- c("2023-01-01", "2023-01-31")

# Retrieve the Sentinel-2 data cube
sentinel2_cube <- get.sentinel2_cube(s_collection, shape, date_range,
                                      aggregation_method = "mean", get.dataset = TRUE)

# Print the resulting data cube details
print(sentinel2_cube)
} # }
```
