#' ToolsRTM: Radiative Transfer Model Tools for Vegetation Remote Sensing
#'
#' A series of tools for simulating spectra with radiative transfer models
#' and Sentinel-2 data.
#'
#' @import ggplot2
#' @importFrom stats aggregate as.formula complete.cases cor dnorm formula IQR
#'   integrate kmeans lm median na.omit na.pass plnorm pnorm prcomp predict
#'   qnorm qpois quantile rlnorm
#' @importFrom utils flush.console read.table setTxtProgressBar stack
#'   txtProgressBar write.csv write.table
#' @importFrom graphics abline axis grid layout legend lines par points rect
#' @importFrom grDevices colorRampPalette rainbow
#' @importFrom pls mvr crossval
#' @importFrom sf st_crs
#' @importFrom dplyr n
#' @keywords internal
"_PACKAGE"
