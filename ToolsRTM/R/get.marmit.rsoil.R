#' Build a canopy-model-ready soil reflectance spectrum from MARMIT
#'
#' High-level wrapper around \code{\link{get.marmit1}}/\code{\link{get.marmit2}}
#' that loads the measured dry-soil reference spectrum and water-optics
#' constants bundled with ToolsRTM, runs the MARMIT model, and resamples the
#' result onto a fixed wavelength grid so it can be dropped straight into
#' \code{rsoil} for \code{\link{foursail}}, \code{\link{foursail2}},
#' \code{\link{inform}}, or into \code{\link{SPART}}'s \code{rsoil} override.
#'
#' Only the \strong{Bablet 2016} soil database (Bablet et al., 2018) is
#' bundled with the package, to keep install size small. Other MARMIT
#' databases (Dupiau 2020, Humper 2015, Lesaignoux 2008, Liu 2002, Lobell
#' 2002, Marcq 2012, Philpot 2014 -- see
#' \url{https://pss-gitlab.math.univ-paris-diderot.fr/marmit/marmit}) can be
#' used by dropping a folder in the same layout (an index CSV
#' \code{<name>/<name>.csv} with columns \code{ID, Refl_file, SMCg, K, a, psi},
#' plus \code{<name>/spectra/<Refl_file>} tab-separated \code{Wvl,R} files)
#' under \code{system.file("extdata", "marmit", "databases", package =
#' "ToolsRTM")} and passing \code{database = "<name>"}.
#'
#' The dry-soil reference for a given \code{id} is the driest spectrum on
#' file for that soil (the row with the smallest \code{SMCg}), matching how
#' the original MARMIT Shiny app selects it. Wavelengths beyond the native
#' range of that spectrum (some databases stop at 2400 or 2490nm, not 2500)
#' are held constant at the last available value -- the same boundary
#' behavior \code{\link{SPART}} uses for its own BSM soil beyond 2400nm.
#'
#' @param database character. Soil database name (a folder under
#'   \code{extdata/marmit/databases/}). Default \code{"Bablet_2016"}, the
#'   only one bundled.
#' @param id integer. Soil ID within the database's index CSV (see the
#'   \code{ID} column). Default \code{1}.
#' @param version character. \code{"marmit1"} (default) or \code{"marmit2"}.
#'   MARMIT-2 additionally accounts for soil particle size/refractive index
#'   (\code{n_i}/\code{k_i}/\code{d_i}) and is generally more accurate for
#'   coarser soils; MARMIT-1 is simpler and matches the original 2018 paper.
#' @param L numeric. Thickness of the surface water layer, cm. Default
#'   \code{0.05}.
#' @param eps numeric (0-1). Fraction of the soil surface that is wet.
#'   Default \code{0.3}.
#' @param n_i,k_i,d_i numeric. MARMIT-2-only soil-particle parameters (real
#'   refractive index, imaginary refractive index, particle volume
#'   fraction). Ignored when \code{version = "marmit1"}. Defaults
#'   (\code{1.53}, \code{0.001}, \code{0.0005}) match the MARMIT Shiny app's
#'   defaults.
#' @param wl.out integer vector. Wavelength grid (nm) the output is
#'   resampled/padded onto. Default \code{400:2500} (matches
#'   \code{\link{foursail}}'s \code{spectrum.all = TRUE} / \code{\link{inform}}
#'   grid); pass \code{400:2400} for Fluspect-leaf-model calls to
#'   \code{\link{foursail}}/\code{\link{foursail2}}.
#'
#' @return A list:
#'   \item{wavelength}{the \code{wl.out} grid.}
#'   \item{rsoil.dry}{dry-soil reflectance on \code{wl.out} (no MARMIT
#'     wetting applied -- the raw measured reference).}
#'   \item{rsoil.wet}{MARMIT-simulated reflectance on \code{wl.out}, ready to
#'     use as \code{rsoil}.}
#'   \item{SMC}{estimated gravimetric soil moisture content (percent),
#'     from \code{\link{sigmoid.soil}}.}
#'   \item{params}{a one-row data.frame recording \code{database}, \code{id},
#'     \code{version}, \code{L}, \code{eps}, \code{n_i}, \code{k_i}, \code{d_i}.}
#'
#' @seealso \code{\link{get.marmit1}}, \code{\link{get.marmit2}}
#' @export
#'
#' @examples
#' \dontrun{
#' soil <- get.marmit.rsoil(database = "Bablet_2016", id = 1, L = 0.05, eps = 0.3)
#' plot(soil$wavelength, soil$rsoil.wet, type = "l")
#' lines(soil$wavelength, soil$rsoil.dry, col = "grey50")
#'
#' # Feed straight into fourSAIL (PROSPECT-D domain, full 400-2500nm)
#' LUT <- as.data.frame(getLUT(inputs = ToolsRTM::inputsPROSAIL, nLUT = 1, setseed = 1))
#' sim <- foursail(inputLUT = LUT, rsoil = soil$rsoil.wet, LeafModel = "PROSPECT-D")
#' }
get.marmit.rsoil <- function(database = "Bablet_2016", id = 1, version = "marmit1",
                              L = 0.05, eps = 0.3, n_i = 1.53, k_i = 0.001, d_i = 0.0005,
                              wl.out = 400:2500) {

  if (!version %in% c("marmit1", "marmit2")) {
    stop("get.marmit.rsoil(): 'version' must be 'marmit1' or 'marmit2', got '", version, "'.")
  }

  db_dir <- system.file("extdata", "marmit", "databases", database, package = "ToolsRTM")
  if (db_dir == "" || !dir.exists(db_dir)) {
    available <- list.dirs(system.file("extdata", "marmit", "databases", package = "ToolsRTM"),
                            full.names = FALSE, recursive = FALSE)
    stop("get.marmit.rsoil(): soil database '", database, "' not found under inst/extdata/marmit/databases. ",
         "Available: ", paste(available, collapse = ", "), ". ",
         "Other MARMIT databases can be added by copying the same folder layout in -- see ?get.marmit.rsoil.")
  }

  index_file <- file.path(db_dir, paste0(database, ".csv"))
  df <- utils::read.csv(index_file)
  df1 <- df[df$ID == id, ]
  if (nrow(df1) == 0) {
    stop("get.marmit.rsoil(): no rows with ID = ", id, " in '", database, "'. ",
         "Available IDs: ", paste(sort(unique(df$ID)), collapse = ", "), ".")
  }

  ## ---- dry-soil reference: the driest (lowest SMCg) spectrum for this ID ----
  dry_row <- df1[which.min(df1$SMCg), ]
  Rd_raw <- utils::read.csv(file.path(db_dir, "spectra", dry_row$Refl_file), sep = "\t")

  wl_native <- seq(max(min(Rd_raw$Wvl), 400), max(Rd_raw$Wvl), by = 1)
  Rd <- stats::approx(Rd_raw$Wvl, Rd_raw$R, xout = wl_native, rule = 2)$y

  ## ---- water-optics constants (Segelstein 1981 / Buiteveld-Kou-Wieliczka), clipped to wl_native ----
  param_dir <- system.file("extdata", "marmit", "parameters", package = "ToolsRTM")
  n_w_raw <- utils::read.csv(file.path(param_dir, "n_segelstein.csv"), sep = "\t")
  alpha_w_raw <- utils::read.csv(file.path(param_dir, "alpha_buikouwie.csv"), sep = "\t")
  n_w <- stats::approx(n_w_raw$Wvl, n_w_raw$n, xout = wl_native, rule = 2)$y
  alpha_w <- stats::approx(alpha_w_raw$Wvl, alpha_w_raw$alpha, xout = wl_native, rule = 2)$y

  K <- as.numeric(unique(dry_row$K))
  a <- as.numeric(unique(dry_row$a))
  psi <- as.numeric(unique(dry_row$psi))

  ## ---- run MARMIT ----
  if (version == "marmit2") {
    Rw <- get.marmit2(n_w, alpha_w, n_i, k_i, Rd, L, eps, d_i, wl_native)
  } else {
    Rw <- get.marmit1(n_w, alpha_w, Rd, L, eps)
  }

  phi <- L * eps
  SMC <- sigmoid.soil(phi, K, a, psi)

  ## ---- resample/pad onto the requested output grid ----
  rsoil.wet <- stats::approx(wl_native, Rw, xout = wl.out, rule = 2)$y
  rsoil.dry <- stats::approx(wl_native, Rd, xout = wl.out, rule = 2)$y

  list(
    wavelength = wl.out,
    rsoil.dry = rsoil.dry,
    rsoil.wet = rsoil.wet,
    SMC = SMC,
    params = data.frame(database = database, id = id, version = version,
                         L = L, eps = eps, n_i = n_i, k_i = k_i, d_i = d_i)
  )
}
