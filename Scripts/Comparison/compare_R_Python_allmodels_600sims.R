## Cross-language robustness comparison: ALL leaf models (PROSPECT-D,
## PROSPECT-PRO, Liberty, Fluspect-B, Fluspect-B-Cx) x ALL canopy models
## (fourSAIL, foursail2, INFORM) supported by both ToolsRTM (R) and
## toolsrtm (Python) -- 15 combinations total. Same shared-LUT approach
## as compare_R_Python_toolsrtm_600sims.R, extended to the full matrix
## following the parameter set proven in vignettes/t04-comparing-models.Rmd.

this_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
root <- normalizePath(file.path(dirname(this_file), "..", ".."))
outdir <- file.path(root, "Scripts/Comparison/_out")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

suppressMessages(devtools::load_all(file.path(root, "ToolsRTM"), quiet = TRUE))

n_samples <- 100
wl_2101 <- 400:2500
wl_2001 <- 400:2400

## Structural + PROSPECT-family part of the LUT: real random sampling from
## the package's own parameter-definition table (same source used by the
## single-combination comparison and by the tutorials).
build_base_lut <- function(seed) {
  as.data.frame(ToolsRTM::getLUT(inputs = ToolsRTM::inputsPROSAIL, nLUT = n_samples, setseed = seed))
}

## Extra per-model columns not covered by inputsPROSAIL, jittered around the
## proven-working central values from t04-comparing-models.Rmd's common_lut.
add_liberty_cols <- function(LUT, seed) {
  set.seed(seed)
  n <- nrow(LUT)
  LUT$cell.d        <- runif(n, 25, 80)
  LUT$inter.c       <- runif(n, 0.02, 0.08)
  LUT$baseline.abs  <- runif(n, 0.0004, 0.0008)
  LUT$leaf.thick    <- runif(n, 1.2, 3)
  LUT$albino.abs    <- runif(n, 0, 3)
  LUT$lign.cell     <- runif(n, 1, 6)
  LUT$Nitrogen      <- runif(n, 0.4, 1.8)
  LUT
}
add_fluspect_cols <- function(LUT, seed) {
  set.seed(seed)
  n <- nrow(LUT)
  LUT$Cs  <- runif(n, 0, 0.03)
  LUT$fqe <- runif(n, 0.005, 0.02)
  LUT$Cx  <- runif(n, 0, 0.2)
  LUT
}
add_foursail2_cols <- function(LUT, seed) {
  set.seed(seed)
  n <- nrow(LUT)
  LUT$fraction_brown <- runif(n, 0, 0.3)
  LUT$diss <- runif(n, 0, 1)
  LUT$Cv   <- runif(n, 0.5, 1)
  LUT$Zeta <- runif(n, 0, 0.5)
  LUT
}
add_inform_cols <- function(LUT, seed) {
  set.seed(seed)
  n <- nrow(LUT)
  LUT$LAIu <- runif(n, 0.2, 1.5)
  LUT$sd   <- runif(n, 300, 900)
  LUT$cd   <- runif(n, 2, 7)
  LUT$h    <- runif(n, 10, 30)
  LUT$skyl <- runif(n, 0.05, 0.3)
  LUT
}

leaf_models   <- c("PROSPECT-D", "PROSPECT-PRO", "Liberty", "Fluspect-B", "Fluspect-B-Cx")
canopy_models <- c("fourSAIL", "foursail2", "INFORM")

run_combo <- function(canopy_model, leaf_model, seed) {
  LUT <- build_base_lut(seed)
  if (leaf_model == "Liberty") LUT <- add_liberty_cols(LUT, seed + 1)
  if (leaf_model %in% c("Fluspect-B", "Fluspect-B-Cx")) LUT <- add_fluspect_cols(LUT, seed + 2)
  if (canopy_model == "foursail2") LUT <- add_foursail2_cols(LUT, seed + 3)
  if (canopy_model == "INFORM") LUT <- add_inform_cols(LUT, seed + 4)

  is_fluspect <- leaf_model %in% c("Fluspect-B", "Fluspect-B-Cx")
  ## Output domain: 2001 pts (400-2400 nm) for the two Fluspect leaf models
  ## (their optipar table only spans that range), 2101 pts otherwise --
  ## true for all 3 canopy models alike (inform.R truncates internally,
  ## see r_understorey[1:2001] / r_sail_inf[1:2001]).
  wl <- if (is_fluspect) wl_2001 else wl_2101
  ## rsoil INPUT domain: fourSAIL/foursail2 take rsoil pre-truncated to the
  ## output domain; INFORM always takes the full 2101-length rsoil (its
  ## understorey reflectance is computed at full resolution and only the
  ## final crown-level result gets truncated for Fluspect leaf models).
  rsoil_in <- if (canopy_model == "INFORM") rep(0.15, length(wl_2101)) else rep(0.15, length(wl))

  refl <- matrix(NA_real_, nrow = n_samples, ncol = length(wl))
  for (i in seq_len(n_samples)) {
    row_i <- LUT[i, , drop = FALSE]
    sim <- tryCatch({
      if (canopy_model == "fourSAIL") {
        suppressMessages(foursail(inputLUT = row_i, rsoil = rsoil_in, LeafModel = leaf_model, spectrum.all = !is_fluspect))$rsot
      } else if (canopy_model == "foursail2") {
        suppressMessages(foursail2(inputLUT = row_i, rsoil = rsoil_in, LeafModel = leaf_model))$rsot
      } else {
        suppressMessages(inform(inputLUT = row_i, rsoil = rsoil_in, LeafModel = leaf_model))
      }
    }, error = function(e) NULL)
    if (!is.null(sim) && length(sim) == length(wl)) refl[i, ] <- sim
  }
  ok <- stats::complete.cases(refl)
  cat("R side,", leaf_model, "+", canopy_model, ":", sum(ok), "/", n_samples,
      "valid, domain =", length(wl), "pts\n")
  list(LUT = LUT, refl = refl, ok = ok, wl = wl)
}

summary_rows <- list()
idx <- 0
for (cm in canopy_models) {
  for (lm in leaf_models) {
    idx <- idx + 1
    res <- run_combo(cm, lm, seed = 1000 + idx)
    slug <- paste0(gsub("[^A-Za-z0-9]", "", cm), "_", gsub("[^A-Za-z0-9]", "", lm))
    write.csv(res$LUT, file.path(outdir, paste0("lut_", slug, ".csv")), row.names = FALSE)
    out <- as.data.frame(res$refl); colnames(out) <- paste0("wl", res$wl)
    out$row <- seq_len(n_samples); out$ok <- res$ok
    write.csv(out, file.path(outdir, paste0("refl_R_", slug, ".csv")), row.names = FALSE)
    summary_rows[[idx]] <- data.frame(canopy_model = cm, leaf_model = lm, slug = slug,
                                       n_valid = sum(res$ok), domain_pts = length(res$wl))
  }
}
write.csv(do.call(rbind, summary_rows), file.path(outdir, "allmodels_R_summary.csv"), row.names = FALSE)
cat("Done, all", idx, "leaf x canopy combinations written to", outdir, "\n")
