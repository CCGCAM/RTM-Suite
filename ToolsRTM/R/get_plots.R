#' Get a plot with simulations
#'
#' @param df  a dataframe
#' @param wave  vector with spectral bands in nm
#'
#' @return a plot with average, mean, percentiles
#' @export
#'
#' @examples
#' # Synthetic example: 20 simulated spectra across 5 bands
#' sim_matrix <- matrix(runif(100, 0.1, 0.5), nrow = 20, ncol = 5)
#' wavelengths <- c(490, 560, 665, 705, 740)
#' get.plots(sim_matrix, wavelengths)
get.plots<- function(df, wave){


  # Convert matrix to data frame
  df <- data.frame(df)
  df$row <- 1:nrow(df)  # Add a row identifier

  # Reshape the data to long format
  df_long <- tidyr::gather(df, key = "band", value = "value", -row)
  # Make 'band' an ordered factor with desired order
  df_long$band <- factor(df_long$band, levels = paste0("X", 1:ncol(df)))
  # Calculate average, 25th percentile, and 50th percentile for each band
  summary_stats <- df_long |>
    dplyr::group_by(band) |>
    dplyr::summarise(
      average = mean(value),
      median = median(value),
      percentile_25 = quantile(value, 0.25),
      percentile_50 = quantile(value, 0.50),
      percentile_75 = quantile(value, 0.75)
    )
  summary_stats$band <-wave

  # Plot using ggplot2
  plot.sim <-ggplot(summary_stats, aes(x = band)) +
    geom_line(aes(y = average), color = "black", linewidth = 0.6) +
    geom_line(aes(y = median),  linetype = "dashed", color = "black", linewidth = 0.6) +
    geom_ribbon(aes(ymin = percentile_25, ymax = percentile_75), linetype = "dashed",fill = "black", alpha = 0.3) +
    labs(
      title = "",
      x = "wavelength (nm)",
      y = "Reflectance"
    ) +
    theme_bw()  +  theme(
      text = element_text(size = 14, face='bold'),  # Increase the text size
      axis.title = element_text(size = 16, face = "bold"),  # Make axis titles bold
      axis.text = element_text(face = "bold"),  # Make axis numbers bold
      plot.title = element_text(face = "bold"),  # Make plot title bold
      plot.subtitle = element_text(size = 14, face='bold')  # Adjust subtitle size
    )

  return(plot.sim)
}
