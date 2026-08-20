#' LUT Inversion Using a Radiative Transfer Model (RTM)
#'
#' This function performs the inversion of a Radiative Transfer Model (RTM) based on observed sensor reflectance values.
#' It compares simulated reflectance values from the RTM with observed values and selects the best matches using different merit functions.
#'
#' @param rfl.sensor A matrix with reflectance values of the observed sensor (e.g., rows representing different observations, columns representing wavelengths).
#' @param rfl.rtm A matrix with reflectance values simulated by the Radiative Transfer (RT) model. The structure of this matrix should be similar to `rfl.sensor`.
#' @param LUT A LUT (Look-Up Table) containing the distribution of biophysical parameters used as input in the RT model (e.g., leaf chlorophyll content, water content).
#' @param wave A vector containing the wavelengths corresponding to the columns of `rfl.sensor` and `rfl.rtm`.
#' @param method The merit function used to evaluate the inversion. Options include:
#' \itemize{
#'   \item \code{'merit-RMSE'}: Root Mean Square Error. The default when no method is provided
#'   \item \code{'merit-NRMSE'}: Normalized RMSE (scaled by the range of observed data).
#'   \item \code{'merit-MAE'}: Mean Absolute Error.
#'   \item \code{'merit-NMB'}: Normalized Mean Bias.
#'   \item \code{'merit-FGE'}: Fractional Gross Error.
#'   \item \code{'merit-DWT'}: RMSE computed on discrete wavelet transform (Haar) coefficients instead of raw reflectance.
#'   \item \code{'merit-1stD'}: RMSE computed on the first derivative of the spectra (finite differences along `wave`) instead of raw reflectance.
#'   \item \code{'merit-custom.metric'}: A custom metric defined by users.

#' }
#' @param nOpt The number of optimal solutions (i.e., the best-matching simulated spectra) to select based on the chosen merit function.
#' @param custom_stat An optional custom statistic function. If provided, this will override the default merit function. The custom function should take two arguments: the simulated and observed values and return a single numeric value (the error or difference metric).
#'

#' @return A list with two elements:
#' \itemize{
#'   \item \code{rfl.b}: A matrix of the best-matching reflectance values selected from `rfl.rtm`.
#'   \item \code{LUT.best}: A data frame containing the corresponding biophysical parameters from the LUT for the best solutions.
#' }
#'
#' @export
#' @importFrom wavelets dwt
#'
#' @examples
#' # Simulated example usage:
#' sensor_data <- matrix(runif(100), nrow = 10, ncol = 10) # Simulated sensor reflectance
#' rtm_data <- matrix(runif(100), nrow = 10, ncol = 10)    # Simulated RTM reflectance
#' lut_table <- data.frame(N = runif(10), Cab = runif(10), Cw = runif(10)) # Simulated LUT
#' wavelengths <- seq(400, 700, length.out = 10)  # Simulated wavelengths
#' result <- get.inversionOpt(sensor_data, rtm_data, lut_table, wavelengths, method = 'merit-RMSE', nOpt = 5)
#' print(result)

get.inversionOpt <- function(rfl.sensor = NULL, rfl.rtm = NULL, LUT = NULL,
                             wave = NULL, method = "merit-RMSE", nOpt = NULL,
                             custom_stat =NULL) {
  # Outputs-------
  nSamples <- dim(LUT)[1]
  Metrics_sim <- list()
  number_id <- list()

  Metrics_ <- list()
  LUT_ <- list()

  Table.metrics <- list()
  Table.metrics.nOpt <- list()
  Table.metrics.mean <- list()
  rfl.best <- list()

  # Error metrics functions
  # Root Mean Square Error
  rmse_f <- function(sim, obs) {
    sqrt(mean((sim - obs)^2, na.rm = TRUE))
  }
  # Normalized Root Mean Square Error
  nrmse_f <- function(sim, obs) {
    range_obs <- max(obs, na.rm = TRUE) - min(obs, na.rm = TRUE)
    if (range_obs == 0) {
      return(NA)  # Avoid division by zero
    }
    rmse_f(sim, obs) / range_obs
  }

  #  Mean Absolute Error
  mae_f <- function(sim, obs) {
    mean(abs(sim - obs), na.rm = TRUE)
  }

  # Normalized Mean Bias
  nmb_f <- function(sim, obs) {
    mean_sim <- mean(sim, na.rm = TRUE)
    mean_obs <- mean(obs, na.rm = TRUE)
    if (mean_obs == 0) {
      return(NA)  # Avoid division by zero
    }
    (mean_sim - mean_obs) / mean_obs
  }

  # Fractional Gross Error
  fge_f <- function(sim, obs) {
    mean(2 * abs(sim - obs) / (sim + obs), na.rm = TRUE)
  }

  # First derivative along `wave` (finite differences, base R only -- no hsdar)
  first_derivative_f <- function(mat, wave) {
    mat <- as.matrix(mat)
    d_wave <- diff(wave)
    d_mat <- t(diff(t(mat)))
    sweep(d_mat, 2, d_wave, "/")
  }

  # Discrete wavelet transform (Haar) coefficients, flattened per row
  dwt_f <- function(mat) {
    mat <- as.matrix(mat)
    t(apply(mat, 1, function(row) {
      as.vector(unlist(wavelets::dwt(row, filter = "haar", fast = TRUE)@W))
    }))
  }

  # Select the appropriate merit function based on the 'method' or custom_stat
  if (!is.null(custom_stat)) {
    merit_function <- custom_stat
    method <- 'merit-custom.metric'
  } else {
    # Check if method is NULL and set a default
    if (is.null(method)) {
      merit_function <- rmse_f  # Default to RMSE if no method is provided
      method <- 'merit-RMSE'
      message("No method provided. Defaulting to RMSE.")
    } else {
      # Define merit function based on method
      merit_function <- switch(method,
                               'merit-RMSE' = rmse_f,
                               'merit-NRMSE' = nrmse_f,
                               'merit-MAE' = mae_f,
                               'merit-NMB' = nmb_f,
                               'merit-FGE' = fge_f,
                               'merit-DWT' = rmse_f,
                               'merit-1stD' = rmse_f,
                               'merit-custom.metric' = custom_stat,
                               stop("Invalid method. Choose from 'merit-RMSE', 'merit-NRMSE', 'merit-MAE', 'merit-NMB', 'merit-FGE', 'merit-DWT', 'merit-1stD', 'merit-custom.metric'")
      )
    }
  }

  # 'merit-DWT' and 'merit-1stD' compare a transformation of the spectra
  # (wavelet coefficients / first derivative) rather than raw reflectance.
  # The transform is applied once upfront; `rfl.best` below still selects
  # from the original (untransformed) `rfl.rtm`.
  rfl.sensor.cmp <- rfl.sensor
  rfl.rtm.cmp <- rfl.rtm
  if (identical(method, 'merit-DWT')) {
    rfl.sensor.cmp <- dwt_f(rfl.sensor)
    rfl.rtm.cmp <- dwt_f(rfl.rtm)
  } else if (identical(method, 'merit-1stD')) {
    rfl.sensor.cmp <- first_derivative_f(rfl.sensor, wave)
    rfl.rtm.cmp <- first_derivative_f(rfl.rtm, wave)
  }

  ###############################################################################
  if (!is.null(method)) {
    message(paste('Merit function using', method, 'is processing'))

    progress_bar <- txtProgressBar(min = 0, max = dim(rfl.sensor)[1], style = 3, char = "=")
    for (i in 1:dim(rfl.sensor)[1]) {
      rfl.sensor.i <- rfl.sensor.cmp[i, ]
      for (j in 1:nSamples) {
        rfl.rtm.i <- rfl.rtm.cmp[j, ]
        # Use the selected merit function
        merit_value <- merit_function(rfl.rtm.i, rfl.sensor.i)
        Metrics_sim[[j]] <- merit_value
        number_id[j] <- j
      }

      Metrics_j <- do.call(rbind, lapply(Metrics_sim, as.data.frame))

      Table.metrics <- cbind(ID_lut = 1:dim(LUT)[1], Metrics_j)
      colnames(Table.metrics) <- c('ID_lut', method)
      Table.metrics <- cbind(Table.metrics, LUT)

      # Order by the chosen merit
      Table.metrics <- Table.metrics[order(Table.metrics[, method]), ]
      best <- Table.metrics$ID_lut[1:nOpt]

      # Select the best (nOpt) solutions
      Table.metrics.nOpt[[i]] <- Table.metrics[1:nOpt, ]
      Table.metrics.mean[[i]] <- round(colMeans(Table.metrics.nOpt[[i]]), 3)

      if (nOpt == 1) {
        rfl.best[[i]] <- rfl.rtm[best[nOpt], ]
      } else {
        rfl.best[[i]] <- colMeans(rfl.rtm[best[1:nOpt], ])
      }

      Metrics_[[i]] <- Metrics_sim

      setTxtProgressBar(progress_bar, value = i)
    }
    close(progress_bar)

    # For the lowest spectra from PROSAIL (nOpt with lowest merit)
    LUT.best <- as.data.frame(do.call(rbind, Table.metrics.mean))
    rfl.b <- do.call(rbind, rfl.best)
    colnames(rfl.b) <- paste('R', wave, sep = '.')
    Table.best <- cbind(ID = c(1:dim(rfl.sensor)[1]), LUT.best, rfl.b)

    return(list(rfl.b, LUT.best))
  }
}
