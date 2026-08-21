#' Spectral global sensitivity analysis (Johnson relative importance per wavelength)
#'
#' Runs a canopy radiative transfer model many times while varying a set of
#' plant/soil traits, then computes the Johnson relative-importance index
#' (\code{\link{get.sobol.indices}}'s \code{I.Johnson_norm} column, via
#' \code{sensitivity::johnson()}) at each wavelength -- i.e. how much each
#' trait relatively contributes to explaining reflectance variance at that
#' wavelength. Produces the data behind a stacked "relative contribution
#' (percent) vs wavelength" plot (traits stacked to 100 percent at every
#' wavelength). Despite the name, this does NOT use \code{get.sobol.indices()}'s
#' own \code{STi} (total Sobol index) column -- that column has a known bug
#' (see \code{ToolsRTM/R/get.sobol.indices.R} and this function's own source
#' comments) and produces a near-uniform, physically-meaningless split across
#' traits; \code{I.Johnson_norm} is the same metric this package's reference
#' sensitivity scripts (\code{TOcheck/sensibilidad/}) already use for this
#' exact figure.
#'
#' @param n.samples integer. Number of RTM simulations to run (also the total
#'   number of rows fed to \code{\link{get.sobol.indices}}, whose own \code{N}
#'   parameter is set to \code{n.samples / 2}, its internal split-half sample
#'   size). Default 1000, matching 500 used for the actual Sobol calculation.
#' @param distribution character. \code{"Uniform"} or \code{"Gaussian"} --
#'   which PDF each varied trait is drawn from (bounds/mean+sd come from
#'   \code{ToolsRTM::inputsPROSAIL}).
#' @param traits character vector. Which trait names (must exist in
#'   \code{ToolsRTM::inputsPROSAIL$variable}) to vary and attribute sensitivity
#'   to. Default: N, Cab, EWT, LMA, LIDFa, LAI -- plus \code{SoilCoef}
#'   (a soil-brightness multiplier, 0.5-1.5, on the flat baseline soil
#'   spectrum) added automatically, matching the classic PROSAIL sensitivity
#'   figure (leaf structure, pigment, water, dry matter, leaf angle, LAI, soil).
#' @param leaf.model,canopy.model character. Passed to \code{\link{simulate_RTM}}.
#' @param rsoil.base numeric. Baseline soil reflectance spectrum (before the
#'   \code{SoilCoef} multiplier); default flat 0.15 across 400-2500nm.
#' @param wl.step integer. Compute Sobol indices every \code{wl.step} nm
#'   instead of at every single nm, for speed (the Sobol calculation itself,
#'   not the RTM simulations, dominates runtime at full 1nm resolution).
#'   Default 5.
#' @param seed integer. Random seed for reproducibility.
#' @param n.cores integer. Number of cores for the parallel simulation loop
#'   (via \code{parallel}/\code{doParallel}, same pattern as
#'   \code{Scripts/Sensibility}). Default \code{parallel::detectCores() - 2}.
#' @param chunk.size integer. Simulations are run in chunks of this size
#'   (default 500) instead of all at once, matching the chunking pattern in
#'   \code{Scripts/1-getSCOPE-v3_withChunck.R}. This keeps memory bounded at
#'   large \code{n.samples} (5000, 20000, ...) and, if \code{save.path} is
#'   given, means a crash partway through only loses the in-progress chunk,
#'   not the whole run.
#' @param save.path character or \code{NULL}. If given, each chunk's raw
#'   simulated reflectance + LUT rows are saved to
#'   \code{file.path(save.path, "chunk_<i>.rds")} as soon as that chunk
#'   finishes (so results survive a crash/interrupt), and are read back from
#'   disk (not re-simulated) if the files already exist -- re-running with
#'   the same \code{save.path} resumes instead of restarting.
#'
#' @return A data.frame with columns \code{wavelength}, \code{trait}, and
#'   \code{STi_pct} (Johnson relative-importance index, normalized to sum to
#'   100% across traits at each wavelength -- the column name is kept as
#'   \code{STi_pct} for backward compatibility even though it is no longer a
#'   Sobol total index, see the source comments) -- long format, ready for
#'   \code{ggplot2::geom_area(position = "stack")}.
#' @export
#'
#' @examples
#' \dontrun{
#' si_uniform <- get.spectral.sensitivity(n.samples = 1000, distribution = "Uniform")
#' si_normal  <- get.spectral.sensitivity(n.samples = 1000, distribution = "Gaussian")
#' }
#' @importFrom foreach %dopar%
get.spectral.sensitivity <- function(n.samples = 1000, distribution = "Uniform",
                                      traits = c("N", "Cab", "EWT", "LMA", "LIDFa", "LAI"),
                                      leaf.model = "PROSPECT-D", canopy.model = "fourSAIL",
                                      rsoil.base = NULL, wl.step = 5, seed = 123,
                                      n.cores = NULL, chunk.size = 500, save.path = NULL) {

  if (is.null(rsoil.base)) rsoil.base <- rep(0.15, 2101)
  wl_all <- 400:2500

  ## ---- 1. Build the sampling table (bounds/mean/sd per trait, chosen PDF) ----
  inputs <- ToolsRTM::inputsPROSAIL
  # use.default = 1 means "vary this trait"; 0 means "hold it at its default"
  # (named the other way round from what it sounds like -- verified empirically).
  inputs$use.default <- 0
  inputs$use.default[inputs$variable %in% traits] <- 1
  inputs$Distribution[inputs$variable %in% traits] <- distribution
  # TypeLidf selects *which* LIDF formula to use (1 = dladgen, 2 = campbell) --
  # it's categorical, not a continuum, so it must never be varied/interpolated.
  inputs$use.default[inputs$variable == "TypeLidf"] <- 0

  # inputsPROSAIL only has Mean_D/Std_D filled in for Cab; fill in the rest
  # (midpoint / (range/4), a standard placeholder sd) so Gaussian sampling
  # works for every trait being varied here.
  needs_gauss <- inputs$variable %in% traits & (inputs$Mean_D == "-" | is.na(inputs$Mean_D))
  inputs$Mean_D[needs_gauss] <- (as.numeric(inputs$lower[needs_gauss]) + as.numeric(inputs$upper[needs_gauss])) / 2
  inputs$Std_D[needs_gauss] <- (as.numeric(inputs$upper[needs_gauss]) - as.numeric(inputs$lower[needs_gauss])) / 4

  LUT <- as.data.frame(ToolsRTM::getLUT(inputs = inputs, nLUT = n.samples, setseed = seed))

  # Soil-brightness multiplier, sampled the same way as the vegetation traits
  set.seed(seed + 1)
  SoilCoef <- if (distribution == "Uniform") {
    stats::runif(n.samples, 0.5, 1.5)
  } else {
    pmax(0.1, stats::rnorm(n.samples, mean = 1, sd = 0.25))
  }
  LUT$SoilCoef <- SoilCoef
  traits <- c(traits, "SoilCoef")

  ## ---- 2. Run the canopy model for every sample, in chunks (parallel) ----
  if (is.null(n.cores)) {
    n.cores <- max(1, parallel::detectCores() - 2)
    if (nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_"))) n.cores <- min(n.cores, 2L)
  }
  cl <- parallel::makeCluster(n.cores)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  if (!is.null(save.path)) dir.create(save.path, showWarnings = FALSE, recursive = TRUE)

  i <- NULL  # avoid "no visible binding" NOTE for the foreach loop variable
  chunk_starts <- seq(1, n.samples, by = chunk.size)
  refl_chunks <- vector("list", length(chunk_starts))

  for (ci in seq_along(chunk_starts)) {
    start_row <- chunk_starts[ci]
    end_row <- min(start_row + chunk.size - 1, n.samples)
    chunk_file <- if (!is.null(save.path)) file.path(save.path, sprintf("chunk_%03d.rds", ci)) else NULL

    if (!is.null(chunk_file) && file.exists(chunk_file)) {
      # Resume: this chunk was already simulated and saved in a previous run.
      refl_chunks[[ci]] <- readRDS(chunk_file)
      message(sprintf("get.spectral.sensitivity: chunk %d/%d (%d-%d) loaded from disk",
                      ci, length(chunk_starts), start_row, end_row))
      next
    }

    rows_i <- start_row:end_row
    refl_list <- foreach::foreach(i = rows_i, .packages = "ToolsRTM") %dopar% {
      rsoil_i <- rsoil.base * SoilCoef[i]
      sim <- tryCatch(
        # suppressMessages(): foursail()/foursail2()/inform() print a "SAIL with
        # <model> is processing" message on every single call, which floods the
        # console across hundreds/thousands of simulations.
        suppressMessages(simulate_RTM(inputLUT = LUT[i, ], rsoil = rsoil_i,
                                      leaf.model = leaf.model, canopy.model = canopy.model)),
        error = function(e) NULL
      )
      if (!is.null(sim) && !is.null(sim$rsot) && length(sim$rsot) == length(wl_all)) sim$rsot else rep(NA_real_, length(wl_all))
    }
    chunk_refl <- do.call(rbind, refl_list)
    if (!is.null(chunk_file)) saveRDS(chunk_refl, chunk_file)
    refl_chunks[[ci]] <- chunk_refl

    message(sprintf("get.spectral.sensitivity: chunk %d/%d (%d-%d) done",
                    ci, length(chunk_starts), start_row, end_row))
  }
  refl <- do.call(rbind, refl_chunks)

  ok <- stats::complete.cases(refl)
  refl <- refl[ok, , drop = FALSE]
  LUT <- LUT[ok, , drop = FALSE]
  n_ok <- nrow(refl)
  message(sprintf("get.spectral.sensitivity: %d/%d simulations usable", n_ok, n.samples))

  N_sobol <- floor(n_ok / 2)

  ## ---- 3. Sobol total index at every wl.step-th wavelength ----
  wl_sel <- seq(1, length(wl_all), by = wl.step)
  results <- vector("list", length(wl_sel))

  for (k in seq_along(wl_sel)) {
    j <- wl_sel[k]
    df_wl <- as.data.frame(LUT[, traits, drop = FALSE])
    df_wl$refl_here <- refl[, j]

    sob <- tryCatch(
      suppressMessages(suppressWarnings(
        ToolsRTM::get.sobol.indices(df_wl, output = "refl_here", N = N_sobol, normalize = TRUE)
      )),
      error = function(e) NULL
    )
    if (is.null(sob)) next

    # get.sobol.indices() returns one row per trait (columns Band/Parameter/Si/STi/...),
    # not a matrix -- match by name, don't index by position.
    #
    # NB (fix): this used to read sob$STi (the "total Sobol index" column).
    # get.sobol.indices()'s own STi formula (ToolsRTM/R/get.sobol.indices.R)
    # has a real bug -- its STi line uses `m.1[, output] * m.1[, output]`
    # (output squared against itself) instead of the input variable at all,
    # so STi comes out as the same near-constant value for every trait in the
    # loop rather than a real per-trait sensitivity estimate. After this
    # function's own normalization to 100%, that produced an artificial
    # equal split across traits (e.g. Cab/EWT/LAI/SoilCoef all landing near
    # 25-33% at every wavelength) with no relationship to the actual physics
    # -- confirmed by comparing against known PROSPECT/leaf-optics absorption
    # features, and against this repo's own reference sensitivity scripts
    # (`TOcheck/sensibilidad/0-Sensibility_Analisis_v2.R` and siblings), which
    # use `sensitivity::johnson()` for this exact spectral-sensitivity figure
    # and never call Sobol at all (their own Sobol/fast99 attempts are
    # commented out, abandoned). get.sobol.indices()'s own I.Johnson_norm
    # column (via sensitivity::johnson(), a *different*, independently
    # correct code path in the same function) is what those reference
    # scripts' figures are actually built from -- switched to that here
    # instead, rather than trying to fix the shared/exported STi formula
    # itself (higher risk -- get.sobol.indices() is exported and its STi
    # column, while wrong, might be relied on elsewhere; not touched here).
    sti <- sob$I.Johnson_norm[match(traits, sob$Parameter)]
    sti[is.na(sti) | sti < 0] <- 0
    sti_pct <- 100 * sti / sum(sti)

    results[[k]] <- data.frame(wavelength = wl_all[j], trait = traits, STi_pct = sti_pct)
  }

  out <- do.call(rbind, results)
  out$distribution <- distribution
  out
}
