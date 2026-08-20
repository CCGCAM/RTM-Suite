## Verification-only script (not a reference export): builds a full, real
## get.RTMf() call (data.rad from getRTMo, data.leafopt from
## get.fluspect_mSCOPE, placeholder-but-real-shaped data.etau/data.etah) to
## confirm the RTMf.R bug fixes (column-recycling + absfs_nl mixup) run
## cleanly and produce sane (finite, non-NaN, physically plausible-sign)
## output. Not wired into scratch_export.R since RTMf isn't ported to
## Python yet -- this is purely an R-side correctness check.
## Repo root, derived from this script's own location (python/scratch/<file>.R)
## rather than hardcoded, so it works regardless of who runs it or from where.
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
root <- normalizePath(file.path(dirname(this_file), "..", ".."))
suppressMessages(devtools::load_all(file.path(root, "ToolsRTM/R"), quiet = TRUE))
suppressMessages(devtools::load_all(file.path(root, "SCOPEinR/R"), quiet = TRUE))

data.spectral <- SCOPEinR::get.spectra.SCOPE(getSpectral = TRUE)

## ---- leaf optics (PROSPECT-D, broadcast to nl layers) + BSM soil, same as scratch_rtmo_export.R ----
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
cat("getRTMo OK. names(data.rad) has tau_dd?", "tau_dd" %in% names(data.rad), "\n")

## ---- leaf fluorescence optics (Mb/Mf) via get.fluspect_mSCOPE, single profile layer ----
mly <- list(nly = 1, pLAI = 1, pCab = Cab, pEWT = EWT, pCar = Car, pLMA = LMA, pCs = 0.0, pN = N)
leafbio_ml <- list(Cx = 0.5, fqe = 0.01, Prot = 0.001, CBC = 0.008, Anth = 1)
data.leafopt <- SCOPEinR::get.fluspect_mSCOPE(mly = mly, spectral = data.spectral, leafbio = leafbio_ml,
                                               soil = list(rs_thermal = 0.06), optipar = NULL, nl = nl,
                                               step = 5, get.plots = FALSE)
cat("get.fluspect_mSCOPE OK. Mb dim:", paste(dim(data.leafopt$Mb), collapse = "x"), "\n")

## ---- placeholder-but-real fluorescence yield (from get.biochemical, one representative call) ----
TDP_bio <- SCOPEinR::define_temp_response_biochem(getTDP = TRUE)
data.meteo_bio <- list(Q = 800, Cs = 400, Temp = 25, eb = 15, Oa = 209, p = 1013)
data.leafbio_bio <- list(Type = "C3", stressfactor = 1, Vcmax25 = 60, BallBerry0 = 0.01, BallBerrySlope = 9,
                          Rdparam = 0.015, Kn0 = 2.48, Knalpha = 2.83, Knbeta = 0.114, TDP = TDP_bio)
data.opts_bio <- data.frame(Value = rep(1, 10), Name = paste0("opt", 1:10))
res_bio <- SCOPEinR::get.biochemical(data.leafbio = data.leafbio_bio, data.meteo = data.meteo_bio,
                                      data.opts = data.opts_bio, fV = 1, get.plots = FALSE)
eta_val <- res_bio$eta
cat("representative eta from get.biochemical:", eta_val, "\n")

data.etau <- array(eta_val, dim = c(nl, 13, 36))
data.etah <- rep(eta_val, nl)

## ---- the actual verification: does the FIXED get.RTMf run cleanly? ----
data.rad_out <- SCOPEinR::get.RTMf(data.spectral, data.rad, data.soil, data.leafopt, data.canopy,
                                    data.gap, data.angles, data.etau, data.etah, get.plots = FALSE)

cat("\nget.RTMf ran without error.\n")
cat("any NA in LoF_?", any(is.na(data.rad_out$LoF_)), "\n")
cat("any NA in EoutF_?", any(is.na(data.rad_out$EoutF_)), "\n")
cat("range LoF_:", range(data.rad_out$LoF_, na.rm = TRUE), "\n")
cat("range EoutF_:", range(data.rad_out$EoutF_, na.rm = TRUE), "\n")
cat("EoutF (integrated):", data.rad_out$EoutF, "\n")
cat("LoutF (integrated):", data.rad_out$LoutF, "\n")
cat("F685:", data.rad_out$F685, " F740:", data.rad_out$F740, "\n")
