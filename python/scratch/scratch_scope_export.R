## Reference export for the get.SCOPE() end-to-end wrapper Python port
## (task 6, scope.py). Calls the REAL, unmodified get.SCOPE() against the
## package's own bundled example LUT_input.csv/setoptions.csv (the exact
## values distilled here) -- against the fixed R source (RTMo.R/RTMf.R/
## RTMz.R bugs already fixed earlier this session).
## Repo root, derived from this script's own location (python/scratch/<file>.R)
## rather than hardcoded, so it works regardless of who runs it or from where.
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
root <- normalizePath(file.path(dirname(this_file), "..", ".."))
outdir <- file.path(root, "python/scopeinpython/tests/refdata")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

suppressMessages(devtools::load_all(file.path(root, "ToolsRTM/R"), quiet = TRUE))
suppressMessages(devtools::load_all(file.path(root, "SCOPEinR/R"), quiet = TRUE))

LUT <- read.csv(file.path(root, "SCOPEinR/inst/input/LUT_input.csv"))
opts <- read.csv(file.path(root, "SCOPEinR/inst/input/setoptions.csv"))
# calc_spectrum_planck (RTMt_planck.R) is NOT ported -- force off so this
# reference case matches what get_scope() actually implements.
opts$Value[opts$Options == "calc_spectrum_planck"] <- 0

set.seed(1)
res <- get.SCOPE(LUT = LUT, n.LUT = 1, options.SCOPE = opts,
                  optipar = SCOPEinR::optipar2017.ProspectD,
                  path.out = file.path(root, "outs"),
                  get.outputs = "ALL", get.plots = FALSE)

sim <- res[[1]]
rad <- sim$data.rad
fluxes <- sim$data.fluxes
canopy <- sim$data.canopy

cat("nlayers:", canopy$nlayers, "\n")
cat("refl[1:3]:", rad$refl[1:3], "\n")
cat("Tcave:", fluxes$Tcave, "  Rntot:", fluxes$Rntot, "\n")

write.csv(data.frame(wave = sim$data.spectral$wlS, refl = rad$refl, Lo_ = rad$Lo_, Eout_ = rad$Eout_),
          file.path(outdir, "ref_scope_spectral.csv"), row.names = FALSE)

write.csv(data.frame(
  nlayers = canopy$nlayers, LAIsunlit = canopy$LAIsunlit, LAIshaded = canopy$LAIshaded,
  Rnctot = fluxes$Rnctot, lEctot = fluxes$lEctot, Hctot = fluxes$Hctot, Actot = fluxes$Actot,
  Tcave = fluxes$Tcave, Rnstot = fluxes$Rnstot, lEstot = fluxes$lEstot, Hstot = fluxes$Hstot,
  Gtot = fluxes$Gtot, Tsave = fluxes$Tsave, Rntot = fluxes$Rntot, lEtot = fluxes$lEtot, Htot = fluxes$Htot,
  Tsu = sim$data.thermal$Tsu, Tsh = sim$data.thermal$Tsh,
  canopyemis = rad$canopyemis, iter.ebal = sim$iter.ebal,
  Pntot_Cab = canopy$Pntot_Cab, Pnsun_Cab = canopy$Pnsun_Cab, Pnsha_Cab = canopy$Pnsha_Cab,
  Ja = canopy$Ja, PNPQ = canopy$PNPQ, fqe = canopy$fqe,
  EoutF = rad$EoutF, F685 = rad$F685, F740 = rad$F740, LoutF = rad$LoutF
), file.path(outdir, "ref_scope_scalars.csv"), row.names = FALSE)

write.csv(data.frame(Tcu = sim$data.thermal$Tcu, Tch = sim$data.thermal$Tch),
          file.path(outdir, "ref_scope_temps.csv"), row.names = FALSE)

cat("scope reference export done.\n")
