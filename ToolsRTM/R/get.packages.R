#' Get the packages
#'
#' This function retrieves the list of R package dependencies for the specified path
#' using `renv::dependencies()`. If `renv` is not installed, it installs `renv` first. THis is designed for getting the packages
#' that the Shiny app needs in this package
#'
#' @param path The directory path to check for dependencies. Defaults to the current directory.
#' @param confirm logical. If \code{TRUE} (the default) and the session is interactive, ask
#'   for confirmation before installing \pkg{renv} or any missing dependency. If \code{FALSE},
#'   never install anything; missing packages are only reported, not installed. Per CRAN policy,
#'   this function never installs packages in a non-interactive session regardless of \code{confirm}.
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
get.packages <- function(path = ".", confirm = TRUE) {

  may_install <- function(what) {
    if (!interactive()) return(FALSE)
    if (!confirm) return(FALSE)
    isTRUE(utils::askYesNo(paste0("Install ", what, "?")))
  }

  # Check if renv is installed; if not, offer to install it
  if (!requireNamespace("renv", quietly = TRUE)) {
    if (may_install("the 'renv' package (required by get.packages())")) {
      utils::install.packages("renv")
    } else {
      stop("get.packages() requires the 'renv' package. Install it with install.packages(\"renv\") and try again.", call. = FALSE)
    }
  }

  # Get the list of dependencies from renv
  deps <- renv::dependencies(path = path)

  # Extract unique package names from the dependencies
  packages <- unique(deps$Package)
  # Define a list of packages to avoid
  packages_to_avoid <- c("hsdar", "gdalUtils")

  packages_to_print <- packages[!packages %in% packages_to_avoid]

  # Filter out packages that are not available
  missing_packages <- setdiff(packages, rownames(utils::installed.packages()))
  missing_packages <- missing_packages[!missing_packages %in% packages_to_avoid]

  if (length(missing_packages) > 0) {
    if (may_install(paste(missing_packages, collapse = ", "))) {
      message("Installing missing packages: ", paste(missing_packages, collapse = ", "))
      utils::install.packages(missing_packages)
    } else {
      message("Missing packages (not installed): ", paste(missing_packages, collapse = ", "))
    }
  } else {
    message("All packages are already installed.")
  }
 return(packages_to_print)
}
