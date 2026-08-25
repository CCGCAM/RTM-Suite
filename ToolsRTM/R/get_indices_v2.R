#' Compute vegetation indices from reflectance spectra (v2, expanded index set)
#'
#' Interpolates reflectance to 1nm resolution and computes a large set of
#' published vegetation indices (VNIR and/or SWIR domain), one row per
#' input spectrum. This is an expanded version of \code{getIndices()} with
#' additional SWIR-domain and red-edge indices.
#'
#' @param data data.frame or matrix of reflectance spectra, one row per sample, with reflectance columns named by wavelength (e.g. "R.400", "R.405", ...).
#' @param pattern.rfl character. Prefix identifying reflectance columns in \code{data}. Default "R.".
#' @param factor numeric. Multiplier applied to reflectance values before computing indices (e.g. to rescale 0-10000 integer reflectance to 0-1). Default 1 (no rescaling).
#' @param spectral.domain character. One of "VNIR" (400-850nm), "SWIR" (800-2550nm), or "VNIR-SWIR" (400-2550nm). Default "VNIR".
#'
#' @return A data.frame of computed indices, one row per input spectrum, with columns that are entirely NA removed.
#' @export
#'
#' @examples
#' # Synthetic reflectance spectrum, VNIR domain
#' wl <- seq(400, 850, by = 5)
#' rfl <- matrix(runif(length(wl), 0.05, 0.5), nrow = 1)
#' colnames(rfl) <- paste0("R.", wl)
#' idx <- get.indices.v2(as.data.frame(rfl), pattern.rfl = "R.", spectral.domain = "VNIR")
get.indices.v2 <-function (data, pattern.rfl = "R.", factor = NULL, spectral.domain = NULL)
{
  if (is.null(factor)) {
    factor <- 1
  }
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("Error: 'data' must be a dataframe or matrix")
  }
  if (is.matrix(data)) {
    data <- as.data.frame(data)
  }
  valid_domains <- c("VNIR", "SWIR", "VNIR-SWIR")
  if (is.null(spectral.domain)) {
    spectral.domain <- "VNIR"
  }
  if (!spectral.domain %in% valid_domains) {
    stop("Error: Please use a correct spectral domain. Options are 'VNIR', 'VNIR-SWIR', and 'SWIR'.")
  }
  print(paste("Estimating indices using", spectral.domain,
              "domain..."))
  reflectance_cols <- grep(paste0("^", pattern.rfl, "[0-9]+(\\.[0-9]+)?$"),
                           names(data), value = TRUE)
  wavelengths <- as.numeric(gsub(paste0(pattern.rfl, "(\\d+\\.?\\d*)"),
                                 "\\1", reflectance_cols))
  if (spectral.domain == "VNIR") {
    selected_cols <- reflectance_cols[wavelengths >= 400 &
                                        wavelengths <= 850]
  }
  else if (spectral.domain == "SWIR") {
    selected_cols <- reflectance_cols[wavelengths >= 800 &
                                        wavelengths <= 2550]
  }
  else if (spectral.domain == "VNIR-SWIR") {
    selected_cols <- reflectance_cols[wavelengths >= 400 &
                                        wavelengths <= 2550]
  }
  if (length(selected_cols) == 0) {
    stop("Error: No reflectance columns fall within the specified spectral domain.")
  }
  min_ <- round(min(wavelengths)) - 10
  max_ <- round(max(wavelengths)) + 10
  range2interpo <- c(min_:max_)
  indices.list = list()
  total = dim(data)[1]
  barProgress <- txtProgressBar(min = 0, max = total, style = 3)
  for (i in c(1:dim(data)[1])) {
    values <- as.numeric(data[i, ]) * factor
    r = signal::interp1(wavelengths, values, range2interpo,
                        extrap = T)
    names(r) <- range2interpo
    indices = c()
    if (spectral.domain %in% "VNIR") {
      indices["NDVI"] <- (r["800"] - r["670"])/(r["800"] +
                                                  r["670"])
      indices["RDVI"] <- (r["800"] - r["670"])/(r["800"] +
                                                  r["670"])^0.5
      indices["SR"] <- r["800"]/r["670"]
      indices["MSR"] <- (r["800"]/r["670"] - 1)/((r["800"]/r["670"])^0.5 +
                                                   1)
      indices["OSAVI"] <- ((1 + 0.16) * (r["800"] - r["670"])/(r["800"] +
                                                                 r["670"] + 0.16))
      indices["MSAVI"] <- 1/2 * (2 * r["800"] + 1 - sqrt(((2 *
                                                             r["800"] + 1)^2) - 8 * (r["800"] - r["670"])))
      indices["MTVI1"] <- 1.2 * (1.2 * (r["800"] - r["550"]) -
                                   2.5 * (r["670"] - r["550"]))
      indices["MTVI2"] <- (1.5 * (1.2 * (r["800"] - r["550"]) -
                                    2.5 * (r["670"] - r["550"])))/sqrt((2 * r["800"] +
                                                                          1)^2 - (6 * r["800"] - 5 * sqrt(r["670"])) -
                                                                         0.5)
      indices["MCARI"] <- ((r["700"] - r["670"]) - 0.2 *
                             (r["700"] - r["550"])) * (r["700"]/r["670"])
      indices["MCARI1"] <- 1.5 * (2.5 * (r["800"] - r["670"]) -
                                    1.3 * (r["800"] - r["550"]))
      indices["MCARI2"] <- (1.5 * (2.5 * (r["800"] - r["670"]) -
                                     1.3 * (r["800"] - r["550"])))/sqrt((2 * r["800"] +
                                                                           1)^2 - (6 * r["800"] - 5 * sqrt(r["670"])) -
                                                                          0.5)
      indices["EVI"] <- 2.5 * (r["800"] - r["670"])/(r["800"] +
                                                       6 * r["670"] - 7.5 * r["400"] + 1)
      indices["LIC1"] <- (r["800"] - r["680"])/(r["800"] +
                                                  r["680"])
      indices["VOG"] <- r["740"]/r["720"]
      indices["VOG2"] <- (r["734"] - r["747"])/(r["715"] +
                                                  r["726"])
      indices["VOG3"] <- (r["734"] - r["747"])/(r["715"] +
                                                  r["720"])
      indices["GM1"] <- r["750"]/r["550"]
      indices["GM2"] <- r["750"]/r["700"]
      indices["TCARI"] <- 3 * ((r["700"] - r["670"]) -
                                 0.2 * (r["700"] - r["550"]) * (r["700"]/r["670"]))
      indices["T.O"] <- as.numeric(indices["TCARI"])/as.numeric(indices["OSAVI"])
      indices["CI"] <- r["750"]/r["710"]
      indices["TVI"] <- 0.5 * (120 * (r["750"] - r["550"]) -
                                 200 * (r["670"] - r["550"]))
      indices["SRPI"] <- r["430"]/r["680"]
      indices["NPQI"] <- (r["415"] - r["435"])/(r["415"] +
                                                  r["435"])
      indices["NPCI"] <- (r["680"] - r["430"])/(r["680"] +
                                                  r["430"])
      indices["CTR1"] <- r["695"]/r["420"]
      indices["CAR"] <- r["515"]/r["570"]
      indices["DCabxc"] <- r["672"]/((r["550"] * (3 * r["708"])))
      indices["DNCabxc"] <- r["860"]/((r["550"] * r["708"]))
      indices["SIPI"] <- (r["800"] - r["445"])/(r["800"] +
                                                  r["680"])
      indices["CRI550"] <- (1/r["510"]) - (1/r["550"])
      indices["CRI700"] <- (1/r["510"]) - (1/r["700"])
      indices["CRI550m"] <- (1/r["515"]) - (1/r["550"])
      indices["CRI700m"] <- (1/r["515"]) - (1/r["700"])
      indices["RCRI550"] <- (1/r["510"]) - (1/r["550"]) *
        r["770"]
      indices["RCRI700"] <- (1/r["510"]) - (1/r["700"]) *
        r["770"]
      indices["PSRI"] <- (r["680"] - r["500"])/r["750"]
      indices["LIC3"] <- r["440"]/r["740"]
      waves_ <- c(w02 = 497, w03 = 560, w04 = 665, w05 = 704,
                  w06 = 740, w07 = 782, w08 = 835, w8A = 865, w11 = 1614,
                  w12 = 2202, w.800 = 800, w.762 = 762)
      indices["CIre"] <- (r["782"]/r["705"]) - 1
      indices["CIrededge"] <- (r["800"]/r["705"]) - 1
      indices["CIgreen"] <- (r["800"]/r["B3"]) - 1
      indices["Chlred.edge"] <- (r["780"]/r["705"])^(-1)
      indices["CVI"] <- (r["800"] * r["665"])/r["665"]^2
      indices["IRECI"] <- (r["780"] - r["665"])/(r["705"]/r["740"])
      indices["REP"] <- 700 + 40 * (((r["665"] + r["780"])/2) -
                                      r["705"])/(r["740"] - r["705"])
      indices["RVI"] <- (r["800"]/r["665"])
      indices["RedEg1"] <- r["705"]/r["665"]
      indices["RedEg2"] <- (r["705"] - r["665"])/(r["705"] +
                                                    r["665"])
      indices["PRI"] <- (r["570"] - r["530"])/(r["570"] +
                                                 r["530"])
      indices["PRI515"] <- (r["515"] - r["530"])/(r["515"] +
                                                    r["530"])
      indices["PRIM1"] <- (r["512"] - r["531"])/(r["512"] +
                                                   r["531"])
      indices["PRIM2"] <- (r["600"] - r["531"])/(r["600"] +
                                                   r["531"])
      indices["PRIM3"] <- (r["670"] - r["531"])/(r["670"] +
                                                   r["531"])
      indices["PRIM4"] <- (r["570"] - r["531"] - r["670"])/(r["570"] +
                                                              r["531"] + r["670"])
      indices["PRIn"] <- as.numeric(indices["PRI"])/(as.numeric(indices["RDVI"]) *
                                                       r["700"]/r["670"])
      indices["PRI_CI"] <- as.numeric(indices["PRI"]) *
        ((r["760"]/r["700"]) - 1)
      indices["B"] <- r["450"]/r["490"]
      indices["G"] <- r["550"]/r["670"]
      indices["R"] <- r["700"]/r["670"]
      indices["BGI1"] <- r["400"]/r["550"]
      indices["BGI2"] <- r["450"]/r["550"]
      indices["BF1"] <- r["400"]/r["410"]
      indices["BF2"] <- r["400"]/r["420"]
      indices["BF3"] <- r["400"]/r["430"]
      indices["BF4"] <- r["400"]/r["440"]
      indices["BF5"] <- r["400"]/r["450"]
      indices["BRI1"] <- r["400"]/r["690"]
      indices["BRI2"] <- r["450"]/r["690"]
      indices["RGI"] <- r["690"]/r["550"]
      indices["RARS"] <- r["746"]/r["513"]
      indices["LIC2"] <- r["440"]/r["690"]
      indices["HI"] <- (r["534"] - r["698"])/(r["534"] +
                                                r["698"]) - 1/2 * r["704"]
      indices["CUR"] <- (r["675"] * r["690"])/(r["683"])^2
      indices["PSSRa"] <- r["800"]/r["680"]
      indices["PSSRb"] <- r["800"]/r["635"]
      indices["PSSRc"] <- r["800"]/r["470"]
      indices["PSNDc"] <- (r["800"] - r["470"])/(r["800"] +
                                                   r["470"])
      waves_ <- c(w02 = 497, w03 = 560, w04 = 665, w05 = 704,
                  w06 = 740, w07 = 782, w08 = 835, w8A = 865, w11 = 1614,
                  w12 = 2202, w.800 = 800, w.762 = 762)
      indices["CR.red.nir.1"] <- r["740"]/(r["800"] + (waves_["w06"] -
                                                         waves_["w.800"]) * (r["782"] - r["800"])/(waves_["w.762"] -
                                                                                                     waves_["w.800"]))
      indices["CR.red.nir.2"] <- r["704"]/(r["782"] + (waves_["w05"] -
                                                         waves_["w07"]) * (r["740"] - r["782"])/(waves_["w06"] -
                                                                                                   waves_["w.762"]))
      indices["CR.red.nir.3"] <- r["704"]/(r["800"] + (waves_["w05"] -
                                                         waves_["w.800"]) * (r["740"] - r["800"])/(waves_["w06"] -
                                                                                                     waves_["w.800"]))
      indices["CR.red.nir.4"] <- r["665"]/(r["800"] + (waves_["w04"] -
                                                         waves_["w.800"]) * (r["704"] - r["800"])/(waves_["w06"] -
                                                                                                     waves_["w.800"]))
      indices["CR.red.nir.5"] <- r["740"]/(r["800"] + (waves_["w06"] -
                                                         waves_["w.800"]) * (r["704"] - r["800"])/(waves_["w05"] -
                                                                                                     waves_["w.800"]))
      indices["CR.red.nir"] <- r["782"]/(r["800"] + (waves_["w07"] -
                                                       waves_["w.800"]) * (r["704"] - r["800"])/(waves_["w05"] -
                                                                                                   waves_["w.800"]))
      indices["CR.red.nir.7"] <- r["762"]/(r["800"] + (waves_["w.762"] -
                                                         waves_["w.800"]) * (r["704"] - r["800"])/(waves_["w05"] -
                                                                                                     waves_["w.800"]))
      indices.list[[i]] = indices
    }
    if (spectral.domain %in% "SWIR" | spectral.domain %in%
        "VNIR-SWIR") {
      indices["GnyLi"] <- ((r["900"] * r["1050"]) - (r["955"] *
                                                       r["1220"]))/((r["900"] * r["1050"]) + (r["955"] *
                                                                                                r["1220"]))
      indices["GnyLi.w850"] <- indices["GnyLi"] <- ((r["850"] *
                                                       r["1050"]) - (r["955"] * r["1220"]))/((r["850"] *
                                                                                                r["1050"]) + (r["955"] * r["1220"]))
      indices["GnyLi.w950"] <- indices["GnyLi"] <- ((r["950"] *
                                                       r["1050"]) - (r["955"] * r["1220"]))/((r["950"] *
                                                                                                r["1050"]) + (r["955"] * r["1220"]))
      indices["CI1"] <- ((r["736"] - r["735"])/1) * (r["990"]/r["720"])
      indices["CI2"] <- ((r["736"] - r["735"])/1) * (r["900"]/r["720"])
      indices["CI2"] <- ((r["736"] - r["735"])/1) * (r["950"]/r["720"])
      indices["MCARI.1510"] <- ((r["700"] - r["1510"]) -
                                  0.2 * (r["700"] - r["550"])) * (r["700"]/r["1510"])
      indices["TCARI.1510"] <- 3 * ((r["700"] - r["1510"]) -
                                      0.2 * (r["700"] - r["550"]) * (r["700"]/r["1510"]))
      indices["OSAVI.1510"] <- ((1 + 0.16) * (r["800"] -
                                                r["1510"])/(r["800"] + r["1510"] + 0.16))
      indices["TCARI/OSAVI.1510"] <- as.numeric(indices["TCARI.1510"])/as.numeric(indices["OSAVI.1510"])
      indices["NRI.1510"] <- (r["1510"] - r["660"])/(r["1510"] +
                                                       r["660"])
      indices["RSI.990.720"] <- r["990"]/r["720"]
      indices["NDNI"] <- (log10(1/r["1510"]) - log10(1/r["1680"]))/(log10(1/r["1510"]) +
                                                                      log10(1/r["1680"]))
      indices["S1080"] <- (r["1080"] - r["660"])/(r["1080"] +
                                                    r["660"])
      indices["S1260"] <- (r["1260"] - r["660"])/(r["1260"] +
                                                    r["660"])
      indices["N1645"] <- (r["1645"] - r["1715"])/(r["1645"] +
                                                     r["1715"])
      indices["N870"] <- (r["870"] - r["1450"])/(r["870"] +
                                                   r["1450"])
      indices["N850.1510"] <- (r["850"] - r["1510"])/(r["850"] +
                                                        r["1510"])
      indices["NN1510"] <- (r["1510"])/(r["850"])
      waves_swir <- c(w02 = 497, w03 = 560, w04 = 665,
                      w05 = 704, w06 = 740, w07 = 782, w08 = 835, w8A = 865,
                      w11 = 1614, w12 = 2202)
      indices["CR.SWIR"] <- r["1614"]/(r["865"] + (waves_swir["w11"] -
                                                     waves_swir["w8A"]) * (r["2202"] - r["865"])/(waves_swir["w12"] -
                                                                                                    waves_swir["w8A"]))
      indices.list[[i]] = indices
    }
    setTxtProgressBar(barProgress, i)
  }
  close(barProgress)
  indices.df <- do.call(rbind, indices.list)
  print(indices.df)
  df.indices <- indices.df
  df.indices = df.indices[, colSums(is.na(df.indices)) != nrow(df.indices)]
  return(indices.df)
}
