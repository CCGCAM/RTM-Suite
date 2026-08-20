## Reference export for the ebal Python port (ebal.py), the final piece
## of the SCOPE thermal energy-balance chain (task 5). Builds a full
## get.ebal() call against the fixed R source (RTMo.R, RTMz.R/RTMf.R
## bugs already fixed earlier this session).
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

## ---- aerodynamic/vertical-profile params needed by ebal specifically ----
## NB: d (zero-plane displacement) and z0m (roughness length) must be
## physically consistent with hc (canopy height, = 2m here, set above) --
## standard rule-of-thumb d=0.65*hc, z0m=0.1*hc -- and z (measurement
## height) must be well above the canopy, else log((h-d)/z0m) etc go
## negative inside get.resistances.
Cd <- 0.3; rwc <- 100; z0m <- 0.1 * hc; d <- 0.65 * hc; leafwidth_ebal <- 0.1
kV <- -0.4
z_meas <- 10; u_wind <- 3
rbs <- 50; rss <- 500

data.canopy <- list(nlayers = nl, litab = litab, lazitab = lazitab, nlazi = 36, nlincl = 13,
                     LAI = LAI, lidf = lidf, xl = xl, hot = hot,
                     Cd = Cd, rwc = rwc, zo = z0m, d = d, hc = hc, leafwidth = leafwidth_ebal, kV = kV)
data.leafbio <- list(rho_thermal = 0.01, tau_thermal = 0.01, emis = 1 - 0.01 - 0.01, kV = kV,
                      Type = "C3", stressfactor = 1, Vcmax25 = 60, BallBerry0 = 0.01, BallBerrySlope = 9,
                      Rdparam = 0.015, Kn0 = 2.48, Knalpha = 2.83, Knbeta = 0.114,
                      TDP = SCOPEinR::define_temp_response_biochem(getTDP = TRUE))
## Build the REAL (Fluspect, non-zero kChlrel/kCarrel) leaf optics FIRST,
## and use it for the getRTMo() call too -- not a separate zeroed-kChlrel
## PROSPECT-D-only object -- so data.rad$Pnh_Cab/Pnu_Cab (net PAR absorbed
## by Cab) come out genuinely non-zero, matching how ebal is actually
## driven (a Q=0 absorbed-PAR edge case makes get.biochemical's Brent
## Ci-solver choke -- confirmed empirically, not chased further since a
## realistic non-zero-Cab leaf is the correct input anyway).
mly <- list(nly = 1, pLAI = 1, pCab = Cab, pEWT = EWT, pCar = Car, pLMA = LMA, pCs = 0.0, pN = N)
leafbio_ml <- list(Cx = 0.5, fqe = 0.01, Prot = 0.001, CBC = 0.008, Anth = 1)
data.leafopt <- SCOPEinR::get.fluspect_mSCOPE(mly = mly, spectral = data.spectral, leafbio = leafbio_ml,
                                               soil = list(rs_thermal = 0.06), optipar = NULL, nl = nl,
                                               step = 5, get.plots = FALSE)
data.leafopt_rtmo <- list(refl = data.leafopt$refl, tran = data.leafopt$tran,
                           kChlrel = data.leafopt$kChlrel, kCarrel = data.leafopt$kCarrel)
optipar <- SCOPEinR::optipar2017.ProspectD
soilpar <- list(BSMBrightness = 0.5, BSMlat = 25, BSMlon = 45)
emp <- list(SMp = 15, SMC = 25, film = 0.015)
spec <- list(GSV = optipar$GSV, Kw = optipar$Kw, nw = optipar$nw)
rsoil_bsm <- as.numeric(SCOPEinR::getBSM(soilpar, spec, emp))
data.soil <- list(rfl.soil = rsoil_bsm, rs_thermal = 0.06, rbs = rbs, rss = rss, GAM = 1000)
data.angles <- list(tts = tts, tto = tto, psi = psi)
Esun_v <- SCOPEinR::Esun_[[1]]
Esky_v <- SCOPEinR::Esky_[[1]]
atmo <- data.frame(wave = data.spectral$wlIrrad, Esky_ = Esky_v, Esun_ = Esun_v)
data.meteo <- list(Ta = 20, Rin = -999, Rli = -999, ea = 15, Ca = 400, p = 1013, z = z_meas, u = u_wind, Oa = 209)
data.opts <- read.csv(file.path(root, "SCOPEinR/inst/input/setoptions.csv"))
data.opts$Value[data.opts$Options == "lite"] <- 1
data.opts$Value[data.opts$Options == "calc_vert_profiles"] <- 0
data.opts$Value[data.opts$Options == "simulation"] <- 0
data.opts$Value[data.opts$Options == "Fluorescence_model"] <- 0
data.opts$Value[data.opts$Options == "MoninObukhov"] <- 1
data.opts$Value[data.opts$Options == "soil_heat_method"] <- 2

outs <- SCOPEinR:::getRTMo(data.spectral, atmo, data.soil, data.leafopt_rtmo, data.canopy,
                            data.leafbio, data.angles, data.meteo, data.opts,
                            canopy.model = "fourSAIL", get.plots = FALSE)
data.rad <- outs$data.rad
data.gap <- outs$data.gap

cat("data.rad$Pnh_Cab:", data.rad$Pnh_Cab, "\n")
cat("data.rad$Pnu_Cab:", data.rad$Pnu_Cab, "\n")
cat("data.rad$Rnhc:", data.rad$Rnhc, "\n")
cat("data.rad$Rnuc:", data.rad$Rnuc, "\n")
cat("any NA in leafbio TDP:", any(is.na(unlist(data.leafbio$TDP))), "\n")
cat("data.canopy$kV:", data.canopy$kV, "\n")

res <- SCOPEinR::get.ebal(data.rad, data.gap, data.meteo, data.soil, data.canopy, data.leafbio,
                           data.leafopt, data.spectral, data.opts, integrate.layer = "layers",
                           k.maxit = 100, get.plots = FALSE)

write.csv(data.frame(Tcu = res$data.thermal$Tcu, Tch = res$data.thermal$Tch),
          file.path(outdir, "ref_ebal_temps.csv"), row.names = FALSE)
write.csv(data.frame(
  counter = res$iter$counter, maxEBercu = res$iter$maxEBercu, maxEBerch = res$iter$maxEBerch,
  maxEBers = res$iter$maxEBers, Tsu = res$data.thermal$Tsu, Tsh = res$data.thermal$Tsh,
  canopyemis = res$data.rad$canopyemis,
  Rnctot = res$data.fluxes$Rnctot, lEctot = res$data.fluxes$lEctot, Hctot = res$data.fluxes$Hctot,
  Actot = res$data.fluxes$Actot, Tcave = res$data.fluxes$Tcave,
  Rnstot = res$data.fluxes$Rnstot, lEstot = res$data.fluxes$lEstot, Hstot = res$data.fluxes$Hstot,
  Gtot = res$data.fluxes$Gtot, Tsave = res$data.fluxes$Tsave,
  Rntot = res$data.fluxes$Rntot, lEtot = res$data.fluxes$lEtot, Htot = res$data.fluxes$Htot
), file.path(outdir, "ref_ebal_scalars.csv"), row.names = FALSE)

cat("ebal reference export done. nl =", nl, " counter =", res$iter$counter, "\n")
cat("maxEBercu/erch/ers:", res$iter$maxEBercu, res$iter$maxEBerch, res$iter$maxEBers, "\n")
