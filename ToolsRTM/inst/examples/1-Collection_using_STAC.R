# Clear the environment
rm(list = ls())

# Load required libraries
library(sf)         # For working with spatial data and shapefiles
library(gdalcubes)  # For working with satellite data cubes
library(rstac)      # For accessing STAC API from Microsoft Planetary Computer
library(dplyr)      # For data manipulation
library(terra)      # For raster data operations
library(leaflet)  # For interactive maps
library(leafem)  # For additional leaflet functionalities
library(htmlwidgets)  # For additional leaflet functionalities
library(ToolsRTM)
link_m <- 'https://planetarycomputer.microsoft.com/api/stac/v1'
link_a <- 'https://earth-search.aws.element84.com/v0'

collections_ <- rstac::stac(link_m,force_version=T) %>%
  collections() %>%
  get_request()
# Extract and print only the IDs of the collections
collection_ids <- sapply(collections_$collections, function(x) x$id)
print(head(collection_ids))
# Filter collections for Sentinel, Landsat, or MODIS
filtered_collections <- collections_$collections[sapply(collections_$collections, function(x) {
  grepl("sentinel| landsat | modis ", x$id, ignore.case = TRUE)
})]

# Extract and print the IDs of the filtered collections
filtered_collection_ids <- sapply(filtered_collections, function(x) x$id)
print(head(filtered_collection_ids))


# Example scenario data (replace this with your real scenario data)
scenario <- sf::st_read("inst/applications/STAC/www/scenario/fungus_infection.gpkg")
scenari_path <- fs::dir_ls("inst/applications/STAC/www/scenario/")
scenario <- sf::st_read(scenari_path[1])
# Set GDAL cubes options
gdalcubes::gdalcubes_options(parallel = 8)



# Input parameters (replace with actual values as needed)
buffer_size <- 300  # Buffer size around the centroid in meters
date_range <- as.Date(c("2023-05-14", "2023-05-16"))  # Date range for data search
cloud_threshold <- 5  # Cloud cover percentage threshold

# Calculate the bounding box using the buffer size around the scenario centroid
bb <- ToolsRTM::get_bounding_box(scenario, buffer_size)

# Print bounding box to verify the result
print(bb)
# Retrieve the list of collections
plot(scenario)
## Example usage
#date_range <- as.Date(c("2023-06-01", "2023-09-30"))
cloud_threshold <- 25

coleccion_names.microsoft <- c('sentinel-2-l2a','landsat-c2-l2','modis-17A2HGF-061','modis-09A1-061','modis-09Q1-061','modis-11A2-061','modis-15A2H-061','modis-15A3H-061')
print(coleccion_names.microsoft)
coleccion_names.aws <- c('sentinel-s2-l2a','sentinel-s2-l2a-cogs')
print(coleccion_names.aws)


satellite_collection <- get.satellite_collection(scenario=scenario, collection=coleccion_names.microsoft[1], 
                                                 cloud_server = 'microsoft', n.limit=1,
                                                 date_range=date_range, cloud_threshold=5, buffer_size = 1500)
head(satellite_collection[[2]])
satellite_collection[[3]]


crs_cube <- "EPSG:3035"  # Use a projected CRS for raster operations

# Define the cube view (spatial extent and resolution)
view <- cube_view(
  srs = crs_cube,
  extent = list(
    t0 = as.character(date_range[[1]]),
    t1 = as.character(date_range[[2]]),
    left = bb["xmin"],
    right = bb["xmax"],
    top = bb["ymax"],
    bottom = bb["ymin"]
  ),
  dx = 20, dy = 20, dt = "P1D",  # 20m resolution, time unit = 1 day
  aggregation = 'first', #"aggregation method
  resampling = 'near'    # Resampling method
)


# Mask to remove clouds and shadows using the SCL band
mask <- image_mask("SCL", values = c(6,7,8, 9,10,11,12))  # Clouds,snow,water bad pixels and shadows

gdalcubes::gdalcubes_options(parallel = parallel::detectCores()-2)
# Create a raster cube for the given collection and view
cube <- raster_cube(satellite_collection[[3]], view,mask)

# Check available bands in the cube
available_bands <- names(cube)  # Adjust based on how you access bands

# Define the preferred order of bands
preferred_order <- c("B02", "B03", "B04", "B05", "B06", "B07", "B08", "B8A","B11","B12",'SCL')
# Reorder the available bands based on preferred order
selected_bands <- preferred_order[preferred_order %in% available_bands]
print(selected_bands)


# Get average raster cube
avg_raster_cube <- get.sentinel2_cube(s_collection=satellite_collection[[3]], shape=scenario,
                                      date_range, aggregation_method = "mean", get.dataset = F)
avg_raster_cube
# Get average raster cube
avg_raster_cube <- get.sentinel2_cube(satellite_collection[[3]], shape=scenario,
                                      date_range, aggregation_method = "first", get.dataset = F)

# Get average raster cube
avg_raster_cube <- get.satellite_cube(satellite_collection, shape=scenario, 
                                      date_range, aggregation_method = "first", get.dataset = F)



plot(avg_raster_cube*1/10000)

# Convert SpatRaster to data frame
df <- as.data.frame(avg_raster_cube, xy = TRUE, na.rm = TRUE)

# Exclude the 'SCL' column, pivot longer, and multiply values by 1/10000
df_long <- df %>%
  dplyr::select(-SCL) %>%  # Exclude the 'SCL' column
  pivot_longer(cols = starts_with("B"), names_to = "band", values_to = "value") %>%
  mutate(value = value * 1/10000)  # Multiply values by 1/10000

# Find the min and max values to set appropriate color scale limits
min_value <- min(df_long$value, na.rm = TRUE)
max_value <- max(df_long$value, na.rm = TRUE)

# Plot with adjusted color scale limits
ggplot(df_long, aes(x = x, y = y, fill = value)) +
  geom_tile() +
  coord_equal() +
  viridis:::scale_fill_viridis(limits = c(min_value, max_value), na.value = "gray90") +  # Set limits and color for NA values
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
# Summary of the raster cube
print(summary(avg_raster_cube))
