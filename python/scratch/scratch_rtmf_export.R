## Reference export for the RTMf Python port (rtmf.py), building on
## scratch_rtmf_verify.R's setup (which only printed sanity-check summary
## stats). Exports BOTH:
##   1. The native-grid (640-850nm, 4nm step, 53 pts) pre-spline
##      quantities (LoF_native, Fhem_native, piLo1..4) -- captured by
##      evaluating get.RTMf's body up to that point via a locally-patched
##      copy of the function that returns them directly, so these can be
##      compared to the Python port at full floating-point precision
##      (nothing spline-related involved yet).
##   2. The normal, fully-interpolated get.RTMf() outputs (LoF_, EoutF_,
##      etc, on the 211-point spectral$wlF grid) -- for a looser-tolerance
##      sanity check, since R's spline (`fmm`) and the Python port's
##      (`not-a-knot`) are a deliberate, documented, small approximation
##      near the two ends of the native grid (see rtmf.py docstring).
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
leafbio_ml <- list(Cx = 0.5, fqe = 0.01, Prot = 0.001, CBC = 0.008, Anth = 1)
data.leafopt <- SCOPEinR::get.fluspect_mSCOPE(mly = mly, spectral = data.spectral, leafbio = leafbio_ml,
                                               soil = list(rs_thermal = 0.06), optipar = NULL, nl = nl,
                                               step = 5, get.plots = FALSE)

TDP_bio <- SCOPEinR::define_temp_response_biochem(getTDP = TRUE)
data.meteo_bio <- list(Q = 800, Cs = 400, Temp = 25, eb = 15, Oa = 209, p = 1013)
data.leafbio_bio <- list(Type = "C3", stressfactor = 1, Vcmax25 = 60, BallBerry0 = 0.01, BallBerrySlope = 9,
                          Rdparam = 0.015, Kn0 = 2.48, Knalpha = 2.83, Knbeta = 0.114, TDP = TDP_bio)
data.opts_bio <- data.frame(Value = rep(1, 10), Name = paste0("opt", 1:10))
res_bio <- SCOPEinR::get.biochemical(data.leafbio = data.leafbio_bio, data.meteo = data.meteo_bio,
                                      data.opts = data.opts_bio, fV = 1, get.plots = FALSE)
eta_val <- res_bio$eta
data.etau <- array(eta_val, dim = c(nl, 13, 36))
data.etah <- rep(eta_val, nl)

## ---- 1. normal get.RTMf() output (fully interpolated, 211-pt grid) ----
data.rad_out <- SCOPEinR::get.RTMf(data.spectral, data.rad, data.soil, data.leafopt, data.canopy,
                                    data.gap, data.angles, data.etau, data.etah, get.plots = FALSE)

write.csv(data.frame(wlF = data.spectral$wlF, LoF_ = data.rad_out$LoF_, EoutF_ = data.rad_out$EoutF_,
                      LoF_sunlit = data.rad_out$LoF_sunlit, LoF_shaded = data.rad_out$LoF_shaded,
                      LoF_scattered = data.rad_out$LoF_scattered, LoF_soil = data.rad_out$LoF_soil,
                      Femliave_ = data.rad_out$Femliave_),
          file.path(outdir, "ref_rtmf_interp.csv"), row.names = FALSE)
write.csv(data.frame(EoutF = data.rad_out$EoutF, LoutF = data.rad_out$LoutF,
                      F685 = data.rad_out$F685, wl685 = data.rad_out$wl685,
                      F740 = data.rad_out$F740, wl740 = data.rad_out$wl740,
                      F684 = data.rad_out$F684, F761 = data.rad_out$F761),
          file.path(outdir, "ref_rtmf_scalars.csv"), row.names = FALSE)

## ---- 2. native-grid pre-spline quantities, via a locally-patched copy of
## get.RTMf's body that returns them directly (same fixed logic, just
## captures piLo1..4/LoF_/Fhem_ before the interp1() calls) ----
body_txt <- deparse(body(SCOPEinR::get.RTMf))
insert_at <- which(grepl("^\\s*Fhem_ <- Fplu_\\[1, \\]\\s*$", body_txt))
stopifnot(length(insert_at) == 1)
body_txt2 <- append(body_txt, "return(list(LoF_native=LoF_, Fhem_native=Fhem_, piLo1=piLo1, piLo2=piLo2, piLo3=piLo3, piLo4=piLo4))", after = insert_at)
get.RTMf.NATIVE <- eval(parse(text = paste(c(
  "function(data.spectral, data.rad, data.soil, data.leafopt, data.canopy, data.gap, data.angles, data.etau, data.etah, get.plots=TRUE) {",
  body_txt2[-c(1, length(body_txt2))], "}"
), collapse = "\n")))
native <- get.RTMf.NATIVE(data.spectral, data.rad, data.soil, data.leafopt, data.canopy,
                           data.gap, data.angles, data.etau, data.etah, get.plots = FALSE)

write.csv(data.frame(wlF_native = seq(640, 850, by = 4), LoF_native = native$LoF_native,
                      Fhem_native = native$Fhem_native, piLo1 = native$piLo1, piLo2 = native$piLo2,
                      piLo3 = native$piLo3, piLo4 = native$piLo4),
          file.path(outdir, "ref_rtmf_native.csv"), row.names = FALSE)

cat("RTMf reference export done. nl =", nl, "\n")
