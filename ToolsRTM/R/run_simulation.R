#' Run a vegetation reflectance/energy-balance simulation, choosing models interactively
#'
#' The single entry point tying ToolsRTM and SCOPEinR together. Two
#' simulation types are supported:
#'
#' \strong{"optical"} -- canopy reflectance only, no atmosphere/energy
#' balance. Runs entirely within ToolsRTM via \code{simulate_RTM()}. Fast,
#' good for building large LUTs for inversion.
#'
#' \strong{"scope"} -- full Soil-Canopy-Observation-Photochemistry-Energy
#' balance simulation (reflectance + fluorescence + radiance + energy
#' fluxes), via SCOPEinR. Requires the SCOPEinR package to be installed --
#' this function checks for it explicitly and gives an informative error
#' (with the real GitLab install command) rather than a cryptic
#' "could not find function" if it's missing, since SCOPEinR is a
#' companion package, not a hard dependency of ToolsRTM.
#'
#' This function does not invent any simulation logic of its own: the
#' "optical" path calls the existing, tested \code{simulate_RTM()}; the
#' "scope" path calls SCOPEinR's own \code{get.SCOPE()} /
#' \code{get.SCOPE.parallel()} using the same argument names confirmed from
#' real teaching scripts using both packages together.
#'
#' @param simulation.type character. "optical" (ToolsRTM only, canopy reflectance) or "scope" (SCOPEinR, full energy balance). Default "optical".
#' @param inputLUT one row (for "optical") or full multi-row LUT (for "scope") of biophysical parameters.
#' @param rsoil soil reflectance spectrum. Required for "optical"; ignored for "scope" (SCOPEinR builds its own from \code{options.SCOPE}).
#' @param leaf.model character. Leaf model name. For "optical": "PROSPECT-PRO", "PROSPECT-D", "Liberty", "Fluspect-B", "Fluspect-B-Cx". For "scope": SCOPEinR's own naming, e.g. "fluspect-CX".
#' @param canopy.model character. Canopy model name. For "optical": "fourSAIL", "foursail2", "INFORM" (the only three real canopy models in ToolsRTM). For "scope": SCOPEinR's own naming, e.g. "fourSAIL".
#' @param optipar optical parameters. For "optical", passed to Fluspect variants if used. For "scope", passed to SCOPEinR (e.g. \code{SCOPEinR::optipar2021.Pro.CX}).
#' @param options.SCOPE SCOPE run-options table (only used when \code{simulation.type = "scope"} -- see SCOPEinR's own documentation for its structure).
#' @param parallel logical. For "scope" only: run in parallel via \code{SCOPEinR::get.SCOPE.parallel()} instead of \code{SCOPEinR::get.SCOPE()}. Default FALSE.
#' @param get.outputs character. For "scope" only: which SCOPE outputs to return, e.g. "ALL". Default "ALL".
#' @param get.plots logical. For "scope" only: whether SCOPEinR should generate its own diagnostic plots. Default FALSE.
#' @param ... additional arguments passed through to the underlying simulate_RTM() / SCOPEinR call.
#'
#' @return For "optical": the output of \code{simulate_RTM()} (reflectance components). For "scope": the output of SCOPEinR's \code{get.SCOPE()}/\code{get.SCOPE.parallel()} (reflectance, fluorescence, radiance, and energy-balance fluxes).
#' @export
#'
#' @examples
#' \dontrun{
#' # Optical-only simulation (ToolsRTM alone, fast, no energy balance)
#' inputs <- ToolsRTM::inputsPROSAIL
#' LUT <- as.data.frame(ToolsRTM::getLUT(inputs = inputs, nLUT = 1, setseed = 1234))
#' sim <- run_simulation(
#'   simulation.type = "optical",
#'   inputLUT = LUT[1, ], rsoil = rsoil,
#'   leaf.model = "PROSPECT-PRO", canopy.model = "fourSAIL"
#' )
#'
#' # Full SCOPE simulation with energy balance (needs SCOPEinR installed)
#' scope_opts <- read.table("Tables/inputs/setoptions.csv", header = TRUE, sep = ",")
#' sim_scope <- run_simulation(
#'   simulation.type = "scope",
#'   inputLUT = LUT, options.SCOPE = scope_opts,
#'   optipar = SCOPEinR::optipar2021.Pro.CX,
#'   leaf.model = "fluspect-CX", canopy.model = "fourSAIL",
#'   parallel = TRUE, get.outputs = "ALL"
#' )
#' }
run_simulation <- function(simulation.type = "optical",
                            inputLUT, rsoil = NULL,
                            leaf.model = "PROSPECT-PRO", canopy.model = "fourSAIL",
                            optipar = NULL, options.SCOPE = NULL,
                            parallel = FALSE, get.outputs = "ALL", get.plots = FALSE,
                            ...) {

  if (!simulation.type %in% c("optical", "scope")) {
    stop("simulation.type must be 'optical' or 'scope' -- got '", simulation.type, "'.")
  }

  if (simulation.type == "optical") {
    if (is.null(rsoil)) {
      stop("simulation.type = 'optical' requires 'rsoil' (a soil reflectance spectrum).")
    }
    return(simulate_RTM(
      inputLUT = inputLUT, rsoil = rsoil,
      leaf.model = leaf.model, canopy.model = canopy.model,
      optipar = optipar, ...
    ))

  } else { # simulation.type == "scope"
    if (!requireNamespace("SCOPEinR", quietly = TRUE)) {
      stop(
        "simulation.type = 'scope' requires the SCOPEinR package, which isn't installed. Install it with:\n",
        "  remotes::install_gitlab(\"caminoccg/scopeinr\", upgrade = \"never\")"
      )
    }
    if (is.null(options.SCOPE)) {
      stop("simulation.type = 'scope' requires 'options.SCOPE' (SCOPE's run-options table).")
    }

    scope_fn <- if (parallel) SCOPEinR::get.SCOPE.parallel else SCOPEinR::get.SCOPE

    return(scope_fn(
      LUT = inputLUT, options.SCOPE = options.SCOPE, optipar = optipar,
      leaf.model = leaf.model, canopy.model = canopy.model,
      get.outputs = get.outputs, get.plots = get.plots, ...
    ))
  }
}
