## Build a self-contained RTMo test case (leaf optics from PROSPECT-D, soil from BSM)
## and dump both the *inputs* (so Python can replicate exactly) and the *outputs* of
## SCOPEinR:::getRTMo, for numerical verification of the Python port.
## Repo root, derived from this script's own location (python/scratch/<file>.R)
## rather than hardcoded, so it works regardless of who runs it or from where.
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
root <- normalizePath(file.path(dirname(this_file), "..", ".."))
outdir <- file.path(root, "python/scratch/_refdata")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

library(ToolsRTM)
library(SCOPEinR)

data.spectral <- SCOPEinR::get.spectra.SCOPE(getSpectral = TRUE)
nwl <- length(data.spectral$wlS)  # 2162

## ---- leaf optics: PROSPECT-D over 400:2400 (2001 pts), replicated over nl layers ----
N <- 1.5; Cab <- 40; Car <- 8; Anth <- 1; Cbrown <- 0.0; EWT <- 0.01; LMA <- 0.009; alpha <- 40
prospect_out <- ToolsRTM:::prospect_DB(N, Cab, Car, Anth, Cbrown, EWT, LMA, alpha)
refl_1 <- prospect_out$refl[1:2001]   # 400:2400
tran_1 <- prospect_out$tran[1:2001]

LAI <- 3
hc <- 2
leafwidth <- 0.1
hot <- leafwidth / hc
LIDFa <- -0.35; LIDFb <- -0.15
tts <- 30; tto <- 0; psi <- 0

nl <- max(2, ceiling(10 * LAI))
xseq <- seq(-1/nl, -1, length.out = nl)
xl <- c(0, xseq)

canopy_angles <- ToolsRTM:::dladgen(LIDFa, LIDFb)  # same as SCOPEinR::leafangles
lidf <- canopy_angles$lidf
litab <- c(seq(5, 75, 10), seq(81, 89, 2))
lazitab <- seq(5, 355, 10)

data.canopy <- list(nlayers = nl, litab = litab, lazitab = lazitab, nlazi = 36,
                     LAI = LAI, lidf = lidf, xl = xl, hot = hot)

data.leafbio <- list(rho_thermal = 0.01, tau_thermal = 0.01)

data.leafopt <- list(
  refl = matrix(rep(refl_1, nl), nrow = nl, byrow = TRUE),
  tran = matrix(rep(tran_1, nl), nrow = nl, byrow = TRUE),
  kChlrel = matrix(0, nrow = nl, ncol = 2001),
  kCarrel = matrix(0, nrow = nl, ncol = 2001)
)

## ---- soil optics: BSM ----
optipar <- SCOPEinR::optipar2017.ProspectD
soilpar <- list(BSMBrightness = 0.5, BSMlat = 25, BSMlon = 45)
emp <- list(SMp = 15, SMC = 25, film = 0.015)
spec <- list(GSV = optipar$GSV, Kw = optipar$Kw, nw = optipar$nw)
rsoil_bsm <- as.numeric(SCOPEinR::getBSM(soilpar, spec, emp))  # 2001 values 400:2400 (optipar is 400:2400? check length)
cat("length rsoil_bsm:", length(rsoil_bsm), "\n")

data.soil <- list(rfl.soil = rsoil_bsm, rs_thermal = 0.06)

data.angles <- list(tts = tts, tto = tto, psi = psi)

## ---- irradiance: default SCOPEinR built-in Esun_/Esky_ ----
Esun_v <- SCOPEinR::Esun_[[1]]
Esky_v <- SCOPEinR::Esky_[[1]]
atmo <- data.frame(wave = data.spectral$wlIrrad, Esky_ = Esky_v, Esun_ = Esun_v)

data.meteo <- list(Ta = 20, Rin = -999, Rli = -999)

## ---- options (data.opts) ----
data.opts <- read.csv("SCOPEinR/inst/input/setoptions.csv")
data.opts$Value[data.opts$Options == "lite"] <- 1
data.opts$Value[data.opts$Options == "calc_vert_profiles"] <- 0
data.opts$Value[data.opts$Options == "simulation"] <- 0

outs <- SCOPEinR:::getRTMo(data.spectral, atmo, data.soil, data.leafopt, data.canopy,
                            data.leafbio, data.angles, data.meteo, data.opts,
                            canopy.model = "fourSAIL", get.plots = FALSE)

data.rad <- outs$data.rad

write.csv(data.frame(wl = data.spectral$wlS,
                      refl = data.rad$refl, rdd = data.rad$rdd, rsd = data.rad$rsd,
                      rdo = data.rad$rdo, rso = data.rad$rso, Lo_ = data.rad$Lo_,
                      Eout_ = data.rad$Eout_, Esun_ = data.rad$Esun_, Esky_ = data.rad$Esky_),
          file.path(outdir, "ref_RTMo_outputs.csv"), row.names = FALSE)

write.csv(data.frame(wl = 400:2400, refl_leaf = refl_1, tran_leaf = tran_1, rsoil = rsoil_bsm),
          file.path(outdir, "ref_RTMo_inputs_spectral.csv"), row.names = FALSE)

scalar_inputs <- data.frame(k = outs$data.gap$k, K = outs$data.gap$K,
                             nl = nl, LAI = LAI, hot = hot, tts = tts, tto = tto, psi = psi,
                             LIDFa = LIDFa, LIDFb = LIDFb)
write.csv(scalar_inputs, file.path(outdir, "ref_RTMo_scalars.csv"), row.names = FALSE)

write.csv(data.frame(Ps = outs$data.gap$Ps, Po = outs$data.gap$Po, Pso = as.numeric(outs$data.gap$Pso), xl = xl),
          file.path(outdir, "ref_RTMo_gap.csv"), row.names = FALSE)

cat("RTMo reference export done. nl =", nl, "\n")
cat("Eouto:", data.rad$Eouto, "\n")
