#' Get the packages
#'
#' This function retrieves the list of R package dependencies for the specified path
#' using `renv::dependencies()`. If `renv` is not installed, it installs `renv` first. THis is designed for getting the packages
#' that the Shiny app needs in this package
#'
#' @param path The directory path to check for dependencies. Defaults to the current directory.
#'
#' @return A character vector of unique package names found in the project dependencies.
#' @export
#'
#' @examples
#' \dontrun{
#' # Get the packages for the current directory
#' get.packages()
#'
#' # Get the packages for a specified path
#' get.packages("path/to/project")
#'
#' }
get.packages <- function(path = ".") {

  # Check if renv is installed; if not, install it
  if (!requireNamespace("renv", quietly = TRUE)) {
    message("renv is not installed. Installing renv...")
    install.packages("renv")
  }

  # Load the renv package
  #library(renv)

  
  
  # Get the list of dependencies from renv
  deps <- renv::dependencies(path = path)

  # Extract unique package names from the dependencies
  packages <- unique(deps$Package)
  # Define a list of packages to avoid
  packages_to_avoid <- c("hsdar", "gdalUtils")

  packages_to_print <- packages[!packages %in% packages_to_avoid]
  
  # Filter out packages that are not available
  missing_packages <- setdiff(packages, rownames(installed.packages()))
  missing_packages <- missing_packages[!missing_packages %in% packages_to_avoid]
  
  if (length(missing_packages) > 0) {
    message("Installing missing packages: ", paste(missing_packages, collapse = ", "))
    install.packages(missing_packages)
  } else {
    message("All packages are already installed.")
  }
 return(packages_to_print)
}
