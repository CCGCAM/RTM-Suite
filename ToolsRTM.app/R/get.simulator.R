#' Launch an interactive RTM Shiny app
#'
#' Launches one of the suite's interactive Shiny applications, bundled under
#' this package's \code{inst/applications/} (not \code{ToolsRTM}'s -- the
#' apps and their UI/data dependencies live in this separate, GitLab-only
#' companion package, see \code{DESCRIPTION}).
#'
#' @param app A character string indicating which simulator to launch.
#'            Acceptable values are:
#'            - `"PROSAIL"`: Launch the PROSAIL simulator.
#'            - `"PROSAIL-BRDF"`: Launch the PROSAIL-BRF simulator.
#'            - `"MARMIT"`: Launch the MARMIT simulator.
#'            - `"SPART"`: Launch the SPART simulator. This option needs the SCOPEinR package
#'            - `"SCOPE"`: Launch the SCOPE simulator. This option needs the SCOPEinR package
#'            - `"getLUT"`: Launch the configuration LUT app.
#'            - `"RTMs"`: Launch the RTM simulator.
#'            - `"STAC"`: Launch the STAC simulator.
#'            - `"Inversion"`: Launch the Inversion module for retrieving plant traits.
#'            - `"default"`: Launch the default simulator (general).
#' @return Launches the Shiny app (invisibly returns \code{NULL}); called for its side effect.
#' @export
#'
#' @examples
#' \dontrun{
#' get.simulator("PROSAIL")
#' get.simulator("SCOPE")
#' get.simulator("MARMIT")
#' get.simulator("getLUT")
#' get.simulator("STAC")
#' }
get.simulator <- function(app = "PROSAIL") {

  valid_apps <- c("PROSAIL", "PROSAIL-BRDF", "MARMIT", "getLUT", "SPART", "SCOPE", "RTMs", "Inversion", "STAC")

  if (!app %in% valid_apps) {
    stop(paste("Invalid app specified. Please choose one of the following options:",
               paste(valid_apps, collapse = ", ")), call. = FALSE)
  }

  if (app %in% c("SCOPE", "SPART") && !"SCOPEinR" %in% rownames(utils::installed.packages())) {
    stop("The SCOPEinR package is required for this app. Please install it from the GitLab repository.", call. = FALSE)
  }

  app_folder <- switch(app,
    "PROSAIL" = "PROSAIL",
    "PROSAIL-BRDF" = "PROSAIL-BRDF",
    "MARMIT" = "MARMIT",
    "SPART" = "SPART",
    "SCOPE" = "SCOPE",
    "getLUT" = "LUTs",
    "STAC" = "STAC",
    "Inversion" = "Inversion",
    "RTMs"
  )
  appDir <- system.file("applications", app_folder, package = "ToolsRTM.app")

  if (appDir == "") {
    stop("Could not find the app directory. Try re-installing `ToolsRTM.app`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}
