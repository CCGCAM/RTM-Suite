#' Simulate canopy reflectance with any supported leaf model + any supported canopy model
#'
#' A single entry point connecting ToolsRTM's leaf models to its canopy
#' models. There are exactly THREE independent canopy models in this
#' package: \code{"fourSAIL"}, \code{"foursail2"}, and \code{"INFORM"}.
#' Everything else you might see referenced internally
#' (\code{foursail.inform}, \code{foursail.inf}, \code{foursail_t_o},
#' \code{foursail_t_s}) are implementation details INFORM uses internally —
#' not standalone models you call directly. An earlier version of this
#' dispatcher incorrectly exposed those as top-level options; that's fixed
#' here.
#'
#' @param inputLUT one row of a LUT (e.g. from \code{getLUT()}).
#' @param rsoil soil reflectance spectrum, as required by the chosen canopy model.
#' @param leaf.model character. One of "PROSPECT-PRO", "PROSPECT-D", "Liberty", "Fluspect-B", "Fluspect-B-Cx" — all 5 work with all 3 canopy models.
#' @param canopy.model character. One of \code{"fourSAIL"}, \code{"foursail2"}, \code{"INFORM"} — the only three real canopy models in this package.
#' @param PROSPECTversion character. Legacy \code{foursail2}-only switch, superseded by \code{leaf.model}; kept for backward compatibility.
#' @param LUT_GB data.frame. Only used when \code{canopy.model = "foursail2"}: reference leaf parameters for the green and brown vegetation fractions (row 1 = green, row 2 = brown). If NULL, foursail2 uses its own built-in defaults.
#' @param optipar optical parameters, passed through when \code{leaf.model} is a Fluspect variant.
#' @param ... additional arguments passed through to the selected canopy model function.
#'
#' @return A list with element \code{rsot} (bi-directional reflectance factor)
#' always present, plus \code{rdot}/\code{rddt}/\code{rsdt} when the selected
#' canopy model provides them (\code{"fourSAIL"}/\code{"foursail2"}: all four;
#' \code{"INFORM"}: only \code{rsot}, since INFORM's forest-level BRF has no
#' separate components -- the others are \code{NULL}).
#' @export
#'
#' @examples
#' \dontrun{
#' inputs <- ToolsRTM::inputsPROSAIL
#' LUT <- as.data.frame(ToolsRTM::getLUT(inputs = inputs, nLUT = 1, setseed = 1234))
#' sim <- simulate_RTM(inputLUT = LUT[1, ], rsoil = rsoil, leaf.model = "PROSPECT-PRO",
#'                      canopy.model = "fourSAIL")
#' }
simulate_RTM <- function(inputLUT, rsoil, leaf.model = 'PROSPECT-PRO', canopy.model = 'fourSAIL',
                          PROSPECTversion = 'PRO', LUT_GB = NULL, optipar = NULL, ...) {

  if (canopy.model == 'fourSAIL') {
    return(foursail(inputLUT = inputLUT, rsoil = rsoil, LeafModel = leaf.model, ...))

  } else if (canopy.model == 'foursail2') {
    return(foursail2(inputLUT = inputLUT, rsoil = rsoil, LUT_GB = LUT_GB,
                      PROSPECTversion = PROSPECTversion, LeafModel = leaf.model, ...))

  } else if (canopy.model == 'INFORM') {
    # inform() itself returns a bare reflectance vector (already the final
    # forest-level BRF, no separate rdot/rddt/rsdt components exist for it),
    # unlike foursail()/foursail2() which return list(rdot=,rsot=,rddt=,rsdt=).
    # Wrapped in the same list shape here (rdot/rddt/rsdt = NULL) so callers
    # can use sim$rsot uniformly regardless of which canopy model was run.
    return(list(rdot = NULL, rsot = inform(inputLUT = inputLUT, rsoil = rsoil, LeafModel = leaf.model, ...),
                rddt = NULL, rsdt = NULL))

  } else {
    stop("Unknown canopy.model: '", canopy.model, "'. This package has exactly three real canopy ",
         "models: 'fourSAIL', 'foursail2', 'INFORM'. (foursail.inform/foursail.inf/foursail_t_o/",
         "foursail_t_s are internal helpers INFORM uses, not standalone options.)")
  }
}
