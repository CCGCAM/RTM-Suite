# Retrieve Satellite Collections

This function retrieves satellite image collections from the Planetary
Computer STAC API based on specified parameters, including bounding box,
date range, and cloud cover threshold.

## Usage

``` r
get.satellite_collection(
  scenario,
  collection,
  cloud_server = "microsoft",
  n.limit = NULL,
  date_range,
  cloud_threshold,
  buffer_size = NULL
)
```

## Arguments

- scenario:

  A character string representing the scenario, which is used to obtain
  the bounding box.

- collection:

  A character string specifying the name of the satellite data
  collection (e.g., "sentinel-2-l2a").

- cloud_server:

  The cloud server to use. Options are "microsoft" or "amazon". Defaults
  to "microsoft".

- n.limit:

  a integer with the maximum number of images during the extractions.

- date_range:

  A length-2 character/Date vector giving the start and end date to
  search (e.g. `c("2023-06-01", "2023-09-30")`).

- cloud_threshold:

  A numeric value representing the maximum allowable cloud cover
  percentage.

- buffer_size:

  A numeric value indicating the buffer size around the bounding box in
  meters (default is 300).

## Value

A filtered image collection based on the specified criteria.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example usage of the get.satellite_collection function
scenario <- "your_scenario_here"  # Replace with your actual scenario
# names ofr every single cloud server
coleccion.microsoft <- c('sentinel-2-l2a', 'landsat-c2-l2', 'modis-17A2HGF-061',
                          'modis-09A1-061', 'modis-09Q1-061', 'modis-11A2-061',
                          'modis-15A2H-061', 'modis-15A3H-061')
coleccion.aws <- c('sentinel-s2-l2a','sentinel-s2-l2a-cogs','landsat-8-l1-c1')
collection <- "sentinel-2-l2a"  # Sentinel-2 Level-2A collection
date_range <- c("2023-06-01", "2023-09-30")  # Date range for January 2023
cloud_threshold <- 20  # Maximum cloud cover of 20%

# Retrieve the satellite collection
satellite_data <- get.satellite_collection(scenario, collection, date_range, cloud_threshold)

# Print the retrieved collection details
print(satellite_data)
} # }
```
