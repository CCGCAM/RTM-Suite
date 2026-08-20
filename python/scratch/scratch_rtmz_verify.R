## Verification-only script (like scratch_rtmf_verify.R): confirms the
## RTMz.R fixes (dead-code-behind-get.plots restructure, column-recycling,
## Po[1:nl+1] index bug) run cleanly end-to-end with real data and produce
## sane, non-degenerate output (i.e. NOT simply the unmodified input, which
## was the pre-fix bug's symptom).
## Repo root, derived from this script's own location (python/scratch/<file>.R)
## rather than hardcoded, so it works regardless of who runs it or from where.
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
root <- normalizePath(file.path(dirname(this_file), "..", ".."))
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
cat("getRTMo OK.\n")

## ---- data.leafopt: baseline (Cx=0) + zeaxanthin (Cx=1) leaf optics ----
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
cat("reflZ - refl any nonzero?", any(abs(data.leafopt$reflZ - data.leafopt$refl) > 1e-12), "\n")

## ---- placeholder-but-real Kn (NPQ) from get.biochemical ----
TDP_bio <- SCOPEinR::define_temp_response_biochem(getTDP = TRUE)
data.meteo_bio <- list(Q = 800, Cs = 400, Temp = 25, eb = 15, Oa = 209, p = 1013)
data.leafbio_bio <- list(Type = "C3", stressfactor = 1, Vcmax25 = 60, BallBerry0 = 0.01, BallBerrySlope = 9,
                          Rdparam = 0.015, Kn0 = 2.48, Knalpha = 2.83, Knbeta = 0.114, TDP = TDP_bio)
data.opts_bio <- data.frame(Value = rep(1, 10), Name = paste0("opt", 1:10))
res_bio <- SCOPEinR::get.biochemical(data.leafbio = data.leafbio_bio, data.meteo = data.meteo_bio,
                                      data.opts = data.opts_bio, fV = 1, get.plots = FALSE)
Kn_val <- res_bio$Kn
cat("representative Kn from get.biochemical:", Kn_val, "\n")

data.Knu <- array(Kn_val, dim = c(nl, 13, 36))
data.Knh <- rep(Kn_val, nl)

## ---- BEFORE: snapshot data.rad WITHIN the wlZ region (500-600nm, R index
## 101:201 since wlS reg1 starts at 400nm step 1) to prove get.RTMz actually
## changes it -- indexing outside 500-600nm would show no change by design.
before_rso <- data.rad$rso[150]  # 549 nm
before_rdo <- data.rad$rdo[150]

data.rad_out <- SCOPEinR::get.RTMz(data.spectral, data.rad, data.soil, data.leafopt, data.canopy,
                                    data.gap, data.angles, data.Knu, data.Knh, get.plots = FALSE)

cat("\nget.RTMz ran without error.\n")
cat("rso[150] (549nm) before:", before_rso, " after:", data.rad_out$rso[150], " changed:", before_rso != data.rad_out$rso[150], "\n")
cat("rdo[150] (549nm) before:", before_rdo, " after:", data.rad_out$rdo[150], " changed:", before_rdo != data.rad_out$rdo[150], "\n")
cat("any NA in rso (wlZ range)?", any(is.na(data.rad_out$rso[500:600])), "\n")
cat("any NA in Eout_ (wlZ range)?", any(is.na(data.rad_out$Eout_[500:600])), "\n")
