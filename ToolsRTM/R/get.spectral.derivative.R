
#' Calculate spectral derivative
#'
#' This function calculates spectral derivatives using either finite approximation or the Savitzky-Golay filter method.
#'
#' @param df A DataFrame containing spectral data. Columns represent wavelengths and spectral values.
#' @param m The order of the derivative. Default is 1.
#' @param method The method to use for derivative calculation. Options are "finApprox" for finite approximation or "sgolay" for Savitzky-Golay filter. Default is "sgolay".
#' @param get.plot a boolena TRUE or FALSE; if TRUE a plot will be done
#' @return A DataFrame with spectral derivatives calculated based on the specified method.
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming you have a DataFrame named 'spectra_df' with columns representing wavelengths and spectral values
#' derived_df <- get.spectral.derivative(spectra_df, m = 2, method = "sgolay")
#' }

get.spectral.derivative <- function(df, m = 1, method = "sgolay",get.plot=T) {

  # Check if the input is a DataFrame
  if (!inherits(df, "data.frame")) {
    stop("Input must be a DataFrame")
  }

  # Perform spectral derivative calculation based on the method
  if (method == "finApprox") {
    # Finite approximation method

    # Assuming 'wavelength' is the column containing wavelengths
    wavelengths <- as.numeric(colnames(df))

    # Calculate the derivatives
    df <- df |>
      dplyr::mutate(dplyr::across(dplyr::everything(), ~ c(diff(.x) / diff(wavelengths), 0)))
    usage_history <- paste(m, ". derivation using finite approximation", sep = "")

  } else if (method == "sgolay") {
    # Savitzky-Golay filter method
    df <- df |>
      dplyr::mutate(dplyr::across(dplyr::everything(), ~ signal::sgolayfilt(.x, m = m)))
    usage_history <- paste(m, ". derivation using Savitzky-Golay filter", sep = "")
  } else {
    stop("Specified method not found")
  }

  if (get.plot == T){

    # Group the files based on the number of rows in the DataFrame
    df$group <- rep(1:ceiling(nrow(df) / 10), each = round(nrow(df)/10,0), length.out = nrow(df))  # Assuming each group contains 10 rows

    # Calculate the average spectral values by band
    avg_spectrum <- df |>
      dplyr::group_by(group) |>
      dplyr::summarise(dplyr::across(dplyr::everything(), mean))   |>
      tidyr::pivot_longer(cols = -group, names_to = "band", values_to = "value") |>
      dplyr::arrange(band, group) |>
      dplyr::mutate(band = factor(band, levels = names(df)))|>
      dplyr::mutate(group = factor(group))


    # Plot the long-format DataFrame
    plot <- ggplot(avg_spectrum, aes(x = band, y = value, color = group)) +
      geom_point() +
      labs(title = "",
           x = "",
           y = "Derivative Plot (avg)",) +
      theme_bw() + theme(legend.position="right",
                           plot.title = element_text(hjust = 0.5, size=10,face="bold"),
                           panel.background = element_rect(fill="white"),
                           plot.background = element_rect(fill = 'white', color = 'white'),
                           legend.key = element_rect(fill = "white", color = "white"),
                           axis.title = element_text(face="bold", size=12),
                           axis.text.y=element_text(hjust = 0.5, size=10,face="bold"),
                           axis.text.x=element_text(angle = 0, vjust = 0.5,hjust = 1, size=10,face="bold"),
                           legend.title=element_blank())
    print(plot)
  }

  # Return the modified DataFrame
  return(df)
}

# Example usage:
# Assuming you have a DataFrame named 'spectra_df' with columns representing wavelengths and spectral values
# derived_df <- spectral_derivative(spectra_df, m = 2, method = "sgolay")
