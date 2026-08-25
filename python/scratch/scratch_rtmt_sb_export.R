## Reference export for the RTMt.sb Python port (rtmt_sb.py, core
## non-obsdir branch only), same RTMo setup pattern as the earlier
## scratch_rtm*_export.R scripts.
## Repo root, derived from this script's own location (python/scratch/<file>.R)
## rather than hardcoded, so it works regardless of who runs it or from where.
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
root <- normalizePath(file.path(dirname(this_file), "..", ".."))
outdir <- file.path(root, "python/scopeinpython/tests/refdata")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

suppressMessages(devtools::load_all(file.path(root, "ToolsRTM/R"), quiet = TRUE))
suppressMessages(devtools::load_all(file.path(root, "SCOPEinR/R"), quiet = TRUE))

data.spectral <- SCOPEinR::get.spectra.SCOPE(getSpectral = TRUE)
N <- 1.5; Cab <- 40; Car <- 8; Anth <- 1; Cbrown <- 0.0; EWT <- 0.01; LMA <- 0.009; alpha <- 40
prospect_out <- ToolsRTM:::prospect_DB(N, Cab, Car, Anth, Cbrown, EWT, LMA, alpha)
refl_1 <- prospect_out$refl[1:2001]
tran_1 <- prospect_out$tran[1:2001]
LAI <- 3; hc <- 2; leafwidth <- 0.1; hot <- leafwidth / hc
LIDFa <- -0.35; LIDFb <- -0.15
tts <- 30; tto <- 0; psi <- 0
nl <- max(2, ceiling(10 * LAI))
xl <- c(0, seq(-1/nl, -1, length.out = nl))
canopy_angles <- ToolsRTM:::dladgen(LIDFa, LIDFb)
lidf <- canopy_angles$lidf
litab <- c(seq(5, 75, 10), seq(81, 89, 2))
lazitab <- seq(5, 355, 10)
data.canopy <- list(nlayers = nl, litab = litab, lazitab = lazitab, nlazi = 36,
                     LAI = LAI, lidf = lidf, xl = xl, hot = hot)
data.leafbio <- list(rho_thermal = 0.01, tau_thermal = 0.01)
data.leafopt_rtmo <- list(
  refl = matrix(rep(refl_1, nl), nrow = nl, byrow = TRUE),
  tran = matrix(rep(tran_1, nl), nrow = nl, byrow = TRUE),
  kChlrel = matrix(0, nrow = nl, ncol = 2001),
  kCarrel = matrix(0, nrow = nl, ncol = 2001)
)
optipar <- SCOPEinR::optipar2017.ProspectD
soilpar <- list(BSMBrightness = 0.5, BSMlat = 25, BSMlon = 45)
emp <- list(SMp = 15, SMC = 25, film = 0.015)
spec <- list(GSV = optipar$GSV, Kw = optipar$Kw, nw = optipar$nw)
rsoil_bsm <- as.numeric(SCOPEinR::getBSM(soilpar, spec, emp))
data.soil <- list(rfl.soil = rsoil_bsm, rs_thermal = 0.06)
data.angles <- list(tts = tts, tto = tto, psi = psi)
Esun_v <- SCOPEinR::Esun_[[1]]
Esky_v <- SCOPEinR::Esky_[[1]]
atmo <- data.frame(wave = data.spectral$wlIrrad, Esky_ = Esky_v, Esun_ = Esun_v)
data.meteo <- list(Ta = 20, Rin = -999, Rli = -999)
data.opts <- read.csv(file.path(root, "SCOPEinR/inst/input/setoptions.csv"))
data.opts$Value[data.opts$Options == "lite"] <- 1
data.opts$Value[data.opts$Options == "calc_vert_profiles"] <- 0
data.opts$Value[data.opts$Options == "simulation"] <- 0
outs <- SCOPEinR:::getRTMo(data.spectral, atmo, data.soil, data.leafopt_rtmo, data.canopy,
                            data.leafbio, data.angles, data.meteo, data.opts,
                            canopy.model = "fourSAIL", get.plots = FALSE)
data.rad <- outs$data.rad
data.gap <- outs$data.gap

## ---- representative leaf/soil temperatures (not from a real ebal
## convergence -- ebal isn't ported yet -- just plausible values to
## exercise get.RTMt.sb's own arithmetic) ----
set.seed(1)
Tcu <- 25 + 3 * sin(seq(0, pi, length.out = nl))  # sunlit leaves warmer, varies with depth
Tch <- rep(22, nl)                                 # shaded leaves cooler, uniform
Tsu <- 30
Tsh <- 24

data.leafopt <- list(refl = data.leafopt_rtmo$refl, tran = data.leafopt_rtmo$tran)
res <- SCOPEinR::get.RTMt.sb(data.rad, data.soil, data.leafbio, data.canopy, data.leafopt,
                              data.gap, Tcu, Tch, Tsu, Tsh, obsdir = FALSE,
                              data.spectral = data.spectral, data.opts = data.opts, get.plots = FALSE)
out_rad <- res$data.rad

write.csv(data.frame(idx = 0:nl, Emint = out_rad$Emint, Eplut = out_rad$Eplut),
          file.path(outdir, "ref_rtmt_sb_profile.csv"), row.names = FALSE)
write.csv(data.frame(layer = 1:nl, Rnuct = out_rad$Rnuct, Rnhct = out_rad$Rnhct),
          file.path(outdir, "ref_rtmt_sb_net.csv"), row.names = FALSE)
write.csv(data.frame(Eoutte = out_rad$Eoutte, Rnust = out_rad$Rnust, Rnhst = out_rad$Rnhst,
                      Tsu = Tsu, Tsh = Tsh),
          file.path(outdir, "ref_rtmt_sb_scalars.csv"), row.names = FALSE)
write.csv(data.frame(Tcu = Tcu, Tch = Tch), file.path(outdir, "ref_rtmt_sb_temps.csv"), row.names = FALSE)

cat("RTMt.sb reference export done. nl =", nl, "\n")
