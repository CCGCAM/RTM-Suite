#' get.RTM - Run Radiative Transfer Model (RTM)
#'
#' This function runs a Radiative Transfer Model (RTM) using specified models for the leaf and canopy, 
#' and allows for optional BRDF correction and plot generation.
#'
#' @param LUT.table A lookup table with input variables for the RTM model (e.g., leaf optical properties, canopy structure).
#' @param leaf.model The leaf model to be used (default is 'fluspect-CX').
#' @param canopy.model The canopy model to be used (default is 'fourSAIL').
#' @param option A string specifying the option or mode (e.g., 'Forward-Mode', 'Inverse-Mode'). Default is 'Forward-Mode'.
#' @param BRDF A boolean to indicate whether Bidirectional Reflectance Distribution Function (BRDF) correction should be applied. Default is TRUE.
#' @param get.plots A boolean to indicate whether to generate plots of the results. Default is TRUE.
#'
#' @return A list with results of the RTM, potentially including reflectance, transmittance, and optional plots.
#' @export
#'
#' @examples
#' \dontrun{
#' # Run RTM with default settings
#' get.RTM(LUT.table)
#'
#' # Run RTM without BRDF correction and plots
#' get.RTM(LUT.table, BRDF = FALSE, get.plots = FALSE)
#' 
#' }
get.RTM<-function(LUT.table,leaf.model='fluspect-CX',canopy.model='fourSAIL',
                    option='Foward-Mode', BRDF=T,
                    get.plots = T) {
  
  
  
  
  ##################################################################################
  ### 0.0 Load optipar (optical leaf properties)
  ##################################################################################
  
  if (missing(LUT.table)){
    stop('please be sure to add a LUT table with main inputs of the RT models')
  } 
  
  ##################################################################################
  ### 0.1 Check for plot or not
  ##################################################################################
  
  if (missing(get.plots)){
    get.plots = FALSE
  }
  
  ##################################################################################
  ### 0.2 Check for leaf model
  ##################################################################################
  
  if (missing(leaf.model)){
 
    stop('please be sure to select a leaf model, e.g., Fluspect-B-Cx, Fluspect-B; Liberty,PROSPECT-D and PROSPECT-PRO')
    
  }
  
  
  ##################################################################################
  ### 0.3 Check for canopy model
  ##################################################################################
  
  if (missing(canopy.model)){
    
    stop('please be sure to select a leaf model, e.g., fourSAIL; fourSAIL2; and INFORM')
    
  }
  
  
  
  
}