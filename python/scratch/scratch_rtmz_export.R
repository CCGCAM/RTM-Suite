## Reference export for the RTMz Python port (rtmz.py), same setup pattern
## as scratch_rtmf_export.R.
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

mly <- list(nly = 1, pLAI = 1, pCab = Cab, pEWT = EWT, pCar = Car, pLMA = LMA, pCs = 0.0, pN = N)
leafbio_ml0 <- list(Cx = 0, fqe = 0.01, Prot = 0.001, CBC = 0.008, Anth = 1)
leafbio_ml1 <- list(Cx = 1, fqe = 0.01, Prot = 0.001, CBC = 0.008, Anth = 1)
data.leafopt <- SCOPEinR::get.fluspect_mSCOPE(mly = mly, spectral = data.spectral, leafbio = leafbio_ml0,
                                               soil = list(rs_thermal = 0.06), optipar = NULL, nl = nl,
                                               step = 5, get.plots = FALSE)
leafoptZ <- SCOPEinR::get.fluspect_mSCOPE(mly = mly, spectral = data.spectral, leafbio = leafbio_ml1,
                                           soil = list(rs_thermal = 0.06), optipar = NULL, nl = nl,
                                           step = 5, get.plots = FALSE)
data.leafopt$reflZ <- leafoptZ$refl
data.leafopt$tranZ <- leafoptZ$tran

TDP_bio <- SCOPEinR::define_temp_response_biochem(getTDP = TRUE)
data.meteo_bio <- list(Q = 800, Cs = 400, Temp = 25, eb = 15, Oa = 209, p = 1013)
data.leafbio_bio <- list(Type = "C3", stressfactor = 1, Vcmax25 = 60, BallBerry0 = 0.01, BallBerrySlope = 9,
                          Rdparam = 0.015, Kn0 = 2.48, Knalpha = 2.83, Knbeta = 0.114, TDP = TDP_bio)
data.opts_bio <- data.frame(Value = rep(1, 10), Name = paste0("opt", 1:10))
res_bio <- SCOPEinR::get.biochemical(data.leafbio = data.leafbio_bio, data.meteo = data.meteo_bio,
                                      data.opts = data.opts_bio, fV = 1, get.plots = FALSE)
Kn_val <- res_bio$Kn
## NB: must be a genuine plain vector (no `dim` attribute) to hit
## get.RTMz's `is.vector(data.Knu)==TRUE` branch -- an array(...,dim=...)
## hits the other branch instead, which (a) has an unverified reshape
## (see rtmz.py's module docstring) and (b) never applies Kn2Cx() at all,
## a further inconsistency not chased down further in this port.
data.Knu <- rep(Kn_val, nl * 13 * 36)
data.Knh <- rep(Kn_val, nl)

## ---- capture the native-grid pre-return-assembly quantities directly:
## get.RTMz doesn't need patching like get.RTMf did, since its return
## value (data.rad) already carries the un-interpolated 500-600nm-band
## corrections at full precision (no spline step in RTMz at all) ----
data.rad_out <- SCOPEinR::get.RTMz(data.spectral, data.rad, data.soil, data.leafopt, data.canopy,
                                    data.gap, data.angles, data.Knu, data.Knh, get.plots = FALSE)

wlZ_native <- 500:600
iwlfi <- match(wlZ_native, data.spectral$wlS)
write.csv(data.frame(
  wl = wlZ_native,
  Lo_before = data.rad$Lo_[iwlfi], Lo_after = data.rad_out$Lo_[iwlfi],
  rso_before = data.rad$rso[iwlfi], rso_after = data.rad_out$rso[iwlfi],
  rdo_before = data.rad$rdo[iwlfi], rdo_after = data.rad_out$rdo[iwlfi],
  Eout_before = data.rad$Eout_[iwlfi], Eout_after = data.rad_out$Eout_[iwlfi]
), file.path(outdir, "ref_rtmz.csv"), row.names = FALSE)

cat("RTMz reference export done. nl =", nl, "\n")
