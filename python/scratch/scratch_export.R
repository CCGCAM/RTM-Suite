## Export reference data + numeric outputs from R for building/verifying the Python port
## Repo root, derived from this script's own location (python/scratch/<file>.R)
## rather than hardcoded, so it works regardless of who runs it or from where.
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
root <- normalizePath(file.path(dirname(this_file), "..", ".."))
outdir <- file.path(root, "python/scratch/_refdata")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## devtools::load_all() (not library()) -- the installed ToolsRTM package can
## lag behind ToolsRTM/R/*.R source fixes (e.g. getIndices() was fixed in
## source without a reinstall), and this script must reflect current source
## behaviour, not a possibly-stale install.
suppressMessages(devtools::load_all(file.path(root, "ToolsRTM/R"), quiet = TRUE))
suppressMessages(devtools::load_all(file.path(root, "SCOPEinR/R"), quiet = TRUE))

## ---- 1. PROSPECT spectral data files ----
write.csv(ToolsRTM::dataSpec_PDB, file.path(outdir, "dataSpec_PDB.csv"), row.names = FALSE)
write.csv(ToolsRTM::dataSpec_PRO, file.path(outdir, "dataSpec_PRO.csv"), row.names = FALSE)

## ---- 2. PROSPECT-D and PROSPECT-PRO reference runs ----
N <- 1.5; Cab <- 40; Car <- 8; Anth <- 1; Cbrown <- 0.2; EWT <- 0.01; LMA <- 0.009; alpha <- 40
db <- ToolsRTM:::prospect_DB(N, Cab, Car, Anth, Cbrown, EWT, LMA, alpha)
write.csv(data.frame(lambda = db$lambda, refl = db$refl, tran = db$tran),
          file.path(outdir, "ref_prospect_DB.csv"), row.names = FALSE)

Prot <- 0.002; CBC <- 0.007
pro <- ToolsRTM::prospect_PRO(N, Cab, Car, Anth, Cbrown, EWT, LMA, alpha, Prot, CBC)
write.csv(data.frame(lambda = pro$lambda, refl = pro$refl, tran = pro$tran),
          file.path(outdir, "ref_prospect_PRO.csv"), row.names = FALSE)

## calctav reference
tav_vals <- data.frame(alpha = c(40, 90, 40, 59),
                        nr = c(1.4, 1.4, 1.5, 1.4))
tav_vals$tav <- mapply(ToolsRTM::calctav, tav_vals$alpha, tav_vals$nr)
write.csv(tav_vals, file.path(outdir, "ref_calctav.csv"), row.names = FALSE)

## ---- 3. fourSAIL reference run (using PROSPECT-PRO leaf model) ----
inputLUT <- data.frame(
  N = N, Cab = Cab, Car = Car, Anth = Anth, Cbrown = Cbrown,
  EWT = EWT, LMA = LMA, alpha = alpha, Prot = Prot, CBC = CBC,
  LIDFa = -0.35, LIDFb = -0.15, TypeLidf = 1, LAI = 3,
  hspot = 0.01, tts = 30, tto = 0, psi = 0
)
rsoil <- rep(0.15, 2101)  # dummy flat soil for 400:2500 (spectrum.all = TRUE => 400:2500, 2101 values)
sail_out <- ToolsRTM::foursail(inputLUT, rsoil, LeafModel = "PROSPECT-PRO", spectrum.all = TRUE)
write.csv(data.frame(lambda = 400:2500, rdot = sail_out$rdot, rsot = sail_out$rsot,
                      rddt = sail_out$rddt, rsdt = sail_out$rsdt),
          file.path(outdir, "ref_foursail.csv"), row.names = FALSE)

## campbell / dladgen reference
camp <- ToolsRTM::campbell(45)
write.csv(data.frame(litab = camp$litab, lidf = camp$lidf), file.path(outdir, "ref_campbell.csv"), row.names = FALSE)
dlad <- ToolsRTM::dladgen(-0.35, -0.15)
write.csv(data.frame(litab = dlad$litab, lidf = dlad$lidf), file.path(outdir, "ref_dladgen.csv"), row.names = FALSE)

## volscatt reference
vs <- ToolsRTM::volscatt(30, 10, 0, 45)
write.csv(data.frame(chi_s = vs[[1]], chi_o = vs[[2]], frho = vs[[3]], ftau = vs[[4]]),
          file.path(outdir, "ref_volscatt.csv"), row.names = FALSE)

cat("ToolsRTM reference export done\n")

## ---- 4. SCOPEinR: constants and BSM data ----
write.csv(SCOPEinR::constants, file.path(outdir, "constants.csv"), row.names = FALSE)

optipar <- SCOPEinR::optipar2017.ProspectD
write.csv(data.frame(wl = optipar$wl, GSV1 = optipar$GSV[,1], GSV2 = optipar$GSV[,2], GSV3 = optipar$GSV[,3],
                      Kw = optipar$Kw, nw = optipar$nw),
          file.path(outdir, "optipar_bsm.csv"), row.names = FALSE)

## ---- 5. BSM reference run ----
soilpar <- list(BSMBrightness = 0.5, BSMlat = 25, BSMlon = 45)
emp <- list(SMp = 15, SMC = 25, film = 0.015)
spec <- list(GSV = optipar$GSV, Kw = optipar$Kw, nw = optipar$nw)
rwet <- SCOPEinR::getBSM(soilpar, spec, emp)
write.csv(data.frame(wl = optipar$wl, rwet = as.numeric(rwet)), file.path(outdir, "ref_BSM.csv"), row.names = FALSE)

## ---- 6. foursail2 reference run (two-layer green/brown canopy, PROSPECT-D) ----
inputLUT_fs2 <- data.frame(
  N = 1.8, Cab = 45, Car = 10, Anth = 1, Cbrown = 0.1,
  EWT = 0.015, LMA = 0.01, alpha = 40, Prot = 0.001, CBC = 0.008,
  LIDFa = 55, LIDFb = 0, TypeLidf = 2,
  LAI = 3, hspot = 0.02, tts = 30, tto = 5, psi = 0,
  fraction_brown = 0.2, diss = 0.5, Cv = 1, Zeta = 0
)
rsoil_fs2 <- rep(0.15, 2101)
fs2_out <- ToolsRTM::foursail2(inputLUT = inputLUT_fs2, rsoil = rsoil_fs2, LeafModel = "PROSPECT-D")
write.csv(data.frame(wl = 400:2500, rdot = fs2_out$rdot, rsot = fs2_out$rsot, rddt = fs2_out$rddt,
                      rsdt = fs2_out$rsdt, alfast = fs2_out$alfast, alfadt = fs2_out$alfadt),
          file.path(outdir, "ref_foursail2.csv"), row.names = FALSE)

## ---- 7. inform reference runs (TypeLidf=1 and TypeLidf=2, PROSPECT-D) ----
make_inform_lut <- function(TypeLidf, LIDFa, LIDFb) {
  data.frame(
    N = 1.8, Cab = 45, Car = 10, Anth = 1, Cbrown = 0.1,
    EWT = 0.015, LMA = 0.01, alpha = 40,
    LIDFa = LIDFa, LIDFb = LIDFb, TypeLidf = TypeLidf,
    LAI = 4, hspot = 0.01, tts = 30, tto = 5, psi = 0,
    LAIu = 0.5, sd = 650, cd = 4.5, h = 20, skyl = 0.1
  )
}
rsoil_inform <- rep(0.15, 2101)
inform_out2 <- ToolsRTM::inform(inputLUT = make_inform_lut(2, 55, 0), rsoil = rsoil_inform, LeafModel = "PROSPECT-D")
write.csv(data.frame(wl = 400:2500, r_forest = inform_out2),
          file.path(outdir, "ref_inform_typelidf2.csv"), row.names = FALSE)
inform_out1 <- ToolsRTM::inform(inputLUT = make_inform_lut(1, -0.35, -0.15), rsoil = rsoil_inform, LeafModel = "PROSPECT-D")
write.csv(data.frame(wl = 400:2500, r_forest = inform_out1),
          file.path(outdir, "ref_inform_typelidf1.csv"), row.names = FALSE)

cat("foursail2 / inform reference export done\n")

## ---- 8. SPART (TOC-only): ToolsRTM's own BSM optipar + a reference TOC run ----
## ToolsRTM::optipar is a separate dataset from SCOPEinR::optipar2017.ProspectD
## (same BSM algorithm, different bundled soil-optics table), so it needs its
## own CSV export -- toolsrtm can't depend on scopeinpython's copy (scopeinpython
## depends on toolsrtm, not the reverse).
write.csv(data.frame(wl = ToolsRTM::optipar$wl, GSV1 = ToolsRTM::optipar$GSV.1,
                      GSV2 = ToolsRTM::optipar$GSV.2, GSV3 = ToolsRTM::optipar$GSV.3,
                      Kw = ToolsRTM::optipar$Kw, nw = ToolsRTM::optipar$nw),
          file.path(outdir, "optipar_spart_bsm.csv"), row.names = FALSE)

inputLUT_spart <- data.frame(
  N = 1.8, Cab = 45, Car = 10, Anth = 1, Cbrown = 0.1,
  EWT = 0.015, LMA = 0.01, alpha = 40, Prot = 0.001, CBC = 0.008,
  LIDFa = 55, LIDFb = 0, TypeLidf = 2,
  LAI = 3, hspot = 0.02, tts = 30, tto = 5, psi = 0
)
data.soil.spart <- list(BSMBrightness = 0.5, BSMlat = 25, BSMlon = 45, rs_thermal = 0.06)
soilemp.spart <- list(SMC = 25, film = 0.015, SMp = 15)
## ToolsRTM::optipar stores GSV as three separate columns (GSV.1/.2/.3), not
## a combined GSV matrix -- SPART() itself reassembles this before calling
## getBSM.toolsRTM (see spart.R); replicate that reassembly here.
GSV.spart <- cbind(ToolsRTM::optipar$GSV.1, ToolsRTM::optipar$GSV.2, ToolsRTM::optipar$GSV.3)
spec.spart <- list(GSV = GSV.spart, Kw = ToolsRTM::optipar$Kw, nw = ToolsRTM::optipar$nw)
rsoil.bsm.spart <- ToolsRTM::getBSM.toolsRTM(soilpar = data.soil.spart, spec = spec.spart, emp = soilemp.spart)
write.csv(data.frame(wl = 400:2400, rsoil = as.numeric(rsoil.bsm.spart)),
          file.path(outdir, "ref_spart_bsm_soil.csv"), row.names = FALSE)

model.sim.spart <- ToolsRTM::foursail(inputLUT = inputLUT_spart[1, ], rsoil = as.numeric(rsoil.bsm.spart),
                                       LeafModel = "PROSPECT-D", spectrum.all = FALSE)
rfl.toc.spart <- ToolsRTM::Compute_BRF(rdot = model.sim.spart$rdot, rsot = model.sim.spart$rsot,
                                        tts = inputLUT_spart[1, "tts"], data.light = ToolsRTM::dataSpec_PDB,
                                        short.waves = TRUE)
write.csv(data.frame(wl = 400:2400, rfl_toc = rfl.toc.spart),
          file.path(outdir, "ref_spart_toc.csv"), row.names = FALSE)

cat("SPART TOC reference export done\n")

## ---- 9. getIndices reference (5 random fourSAIL-simulated spectra) ----
set.seed(42)
wl_idx <- 400:2500
n_idx <- 5
mat_idx <- matrix(NA, nrow = n_idx, ncol = length(wl_idx))
for (i in 1:n_idx) {
  lut_i <- data.frame(
    N = runif(1, 1.2, 2.5), Cab = runif(1, 10, 70), Car = runif(1, 5, 15),
    Anth = runif(1, 0, 2), Cbrown = runif(1, 0, 0.3),
    EWT = runif(1, 0.005, 0.03), LMA = runif(1, 0.003, 0.015), alpha = 40,
    LIDFa = runif(1, -0.5, 0.5), LIDFb = runif(1, -0.3, 0.3), TypeLidf = 1,
    LAI = runif(1, 0.5, 6), hspot = 0.01, tts = 30, tto = 0, psi = 0
  )
  sail_i <- ToolsRTM::foursail(inputLUT = lut_i, rsoil = rep(0.15, 2101), LeafModel = "PROSPECT-D", spectrum.all = TRUE)
  mat_idx[i, ] <- sail_i$rsot
}
colnames(mat_idx) <- paste0("R.", wl_idx)
df_idx <- as.data.frame(mat_idx)
write.csv(df_idx, file.path(outdir, "indices_input_reflectance.csv"), row.names = FALSE)

for (dom in c("VNIR", "SWIR", "VNIR-SWIR")) {
  out_idx <- ToolsRTM::getIndices(df_idx, pattern.rfl = "R.", spectral.domain = dom)
  idx_only <- out_idx[, !(names(out_idx) %in% names(df_idx)), drop = FALSE]
  fname_idx <- paste0("ref_indices_", gsub("-", "_", dom), ".csv")
  write.csv(idx_only, file.path(outdir, fname_idx), row.names = FALSE)
}
cat("getIndices reference export done\n")

## ---- 10. Fluspect-B / Fluspect-B-Cx reference (needed by inform()'s
## Fluspect leaf-model branches and by RTMf) ----
## NB: `optipar` was reassigned to SCOPEinR::optipar2017.ProspectD in
## section 5 above -- must use ToolsRTM::optipar explicitly here, not the
## shadowed bare name (same reasoning as section 8's own comment).
write.csv(ToolsRTM::optipar[, c("wl","nr","Kab","Kca","Ks","Kw","Kdm","phiI","phiII","KcaV","KcaZ","Kant","Kp","Kcbc")],
          file.path(outdir, "optipar_fluspect.csv"), row.names = FALSE)

lut_flus_d <- data.frame(Cab = 40, Car = 10, EWT = 0.015, LMA = 0.01, Cs = 0.1, N = 1.8, fqe = 0.01, Cx = 0.5)
res_flus_d <- ToolsRTM::getFluspect.B(inputsLeaf = lut_flus_d, inputsOptipar = ToolsRTM::optipar, version = 'D')
write.csv(data.frame(wl = res_flus_d$lambda, refl = res_flus_d$refl, tran = res_flus_d$tran,
                      kChlrel = res_flus_d$kChlrel),
          file.path(outdir, "ref_fluspect_b_d_refltran.csv"), row.names = FALSE)
write.csv(res_flus_d$MbI,  file.path(outdir, "ref_fluspect_b_d_MbI.csv"),  row.names = FALSE)
write.csv(res_flus_d$MbII, file.path(outdir, "ref_fluspect_b_d_MbII.csv"), row.names = FALSE)
write.csv(res_flus_d$MfI,  file.path(outdir, "ref_fluspect_b_d_MfI.csv"),  row.names = FALSE)
write.csv(res_flus_d$MfII, file.path(outdir, "ref_fluspect_b_d_MfII.csv"), row.names = FALSE)

## NB: getFluspect.Cx() has a real indexing bug in its fluorescence-matrix
## computation (Iwlf <- intersect(wlp, wlf) used positionally instead of
## via which()/match()) -- this reference run captures R's actual (buggy)
## output on purpose, so the Python port can be verified to match R
## exactly, bug included. See python/toolsrtm/src/toolsrtm/fluspect.py.
lut_flus_cx <- data.frame(Cab = 40, Car = 10, EWT = 0.015, LMA = 0.01, Cs = 0.1, N = 1.8, fqe = 0.01, Cx = 0.5,
                           Prot = 0.001, CBC = 0.008, Anth = 1)
res_flus_cx <- ToolsRTM::getFluspect.Cx(inputsLeaf = lut_flus_cx, inputsOptipar = ToolsRTM::optipar, version = 'Cx')
write.csv(data.frame(wl = res_flus_cx$lambda, refl = res_flus_cx$refl, tran = res_flus_cx$tran,
                      kChlrel = res_flus_cx$kChlrel, kCarrel = res_flus_cx$kCarrel),
          file.path(outdir, "ref_fluspect_cx_refltran.csv"), row.names = FALSE)
write.csv(res_flus_cx$Mb, file.path(outdir, "ref_fluspect_cx_Mb.csv"), row.names = FALSE)
write.csv(res_flus_cx$Mf, file.path(outdir, "ref_fluspect_cx_Mf.csv"), row.names = FALSE)

cat("Fluspect reference export done\n")

## ---- 11. LIBERTY leaf model reference ----
write.csv(ToolsRTM::dataspec.liberty, file.path(outdir, "dataspec_liberty.csv"), row.names = FALSE)

lut_liberty <- data.frame(cell.d = 40, inter.c = 0.045, baseline.abs = 0.0004, leaf.thick = 1.6,
                           albino.abs = 2, Cab = 40, EWT = 0.015, lign.cell = 4, Nitrogen = 1.2)
res_liberty <- ToolsRTM::liberty(inputLUT = lut_liberty)
write.csv(data.frame(wl = res_liberty$lambda, refl = res_liberty$refl, tran = res_liberty$tran, RR = res_liberty$RR),
          file.path(outdir, "ref_liberty.csv"), row.names = FALSE)

cat("LIBERTY reference export done\n")

## ---- 12. Leaf-model wiring reference (Liberty/Fluspect-B/-Cx through
## foursail/foursail2/inform's LeafModel dispatch, end to end) ----
rsoil_2101 <- rep(0.15, 2101)

lut_wire_lib <- data.frame(cell.d = 40, inter.c = 0.045, baseline.abs = 0.0004, leaf.thick = 1.6,
                            albino.abs = 2, Cab = 40, EWT = 0.015, lign.cell = 4, Nitrogen = 1.2,
                            LIDFa = 45, LIDFb = 0, TypeLidf = 2, LAI = 3, hspot = 0.01, tts = 30, tto = 5, psi = 0)
res_wire_lib <- ToolsRTM::foursail(inputLUT = lut_wire_lib, rsoil = rsoil_2101, LeafModel = "Liberty", spectrum.all = TRUE)
write.csv(data.frame(wl = 400:2500, rdot = res_wire_lib$rdot, rsot = res_wire_lib$rsot,
                      rddt = res_wire_lib$rddt, rsdt = res_wire_lib$rsdt),
          file.path(outdir, "ref_foursail_liberty.csv"), row.names = FALSE)

lut_wire_fb <- data.frame(Cab = 40, Car = 10, EWT = 0.015, LMA = 0.01, Cs = 0.1, N = 1.8, fqe = 0.01, Cx = 0.5,
                           LIDFa = 45, LIDFb = 0, TypeLidf = 2, LAI = 3, hspot = 0.01, tts = 30, tto = 5, psi = 0)
res_wire_fb <- ToolsRTM::foursail(inputLUT = lut_wire_fb, rsoil = rsoil_2101, LeafModel = "Fluspect-B", spectrum.all = TRUE)
write.csv(data.frame(wl = 400:2400, rdot = res_wire_fb$rdot, rsot = res_wire_fb$rsot,
                      rddt = res_wire_fb$rddt, rsdt = res_wire_fb$rsdt),
          file.path(outdir, "ref_foursail_fluspectb.csv"), row.names = FALSE)

lut_wire_fcx <- data.frame(Cab = 40, Car = 10, EWT = 0.015, LMA = 0.01, Cs = 0.1, N = 1.8, fqe = 0.01, Cx = 0.5,
                            Prot = 0.001, CBC = 0.008, Anth = 1,
                            LIDFa = 45, LIDFb = 0, TypeLidf = 2, LAI = 3, hspot = 0.01, tts = 30, tto = 5, psi = 0)
res_wire_fcx <- ToolsRTM::foursail(inputLUT = lut_wire_fcx, rsoil = rsoil_2101, LeafModel = "Fluspect-B-Cx", spectrum.all = TRUE)
write.csv(data.frame(wl = 400:2400, rdot = res_wire_fcx$rdot, rsot = res_wire_fcx$rsot,
                      rddt = res_wire_fcx$rddt, rsdt = res_wire_fcx$rsdt),
          file.path(outdir, "ref_foursail_fluspectcx.csv"), row.names = FALSE)

lut_wire_lib2 <- data.frame(cell.d = 40, inter.c = 0.045, baseline.abs = 0.0004, leaf.thick = 1.6,
                             albino.abs = 2, Cab = 40, EWT = 0.015, lign.cell = 4, Nitrogen = 1.2,
                             LIDFa = 45, LIDFb = 0, TypeLidf = 2, LAI = 3, hspot = 0.02, tts = 30, tto = 5, psi = 0,
                             fraction_brown = 0.2, diss = 0.5, Cv = 1, Zeta = 0)
res_wire_lib2 <- ToolsRTM::foursail2(inputLUT = lut_wire_lib2, rsoil = rsoil_2101, LeafModel = "Liberty")
write.csv(data.frame(wl = 400:2500, rdot = res_wire_lib2$rdot, rsot = res_wire_lib2$rsot, rddt = res_wire_lib2$rddt,
                      rsdt = res_wire_lib2$rsdt, alfast = res_wire_lib2$alfast, alfadt = res_wire_lib2$alfadt),
          file.path(outdir, "ref_foursail2_liberty.csv"), row.names = FALSE)

inform_common <- list(LIDFa = 55, LIDFb = 0, TypeLidf = 2, LAI = 4, hspot = 0.01, tts = 30, tto = 5, psi = 0,
                       LAIu = 0.5, sd = 650, cd = 4.5, h = 20, skyl = 0.1)

lut_wire_inf_lib <- data.frame(c(list(cell.d = 40, inter.c = 0.045, baseline.abs = 0.0004, leaf.thick = 1.6,
                                       albino.abs = 2, Cab = 40, EWT = 0.015, lign.cell = 4, Nitrogen = 1.2,
                                       alpha = 40), inform_common))
res_wire_inf_lib <- ToolsRTM::inform(inputLUT = lut_wire_inf_lib, rsoil = rsoil_2101, LeafModel = "Liberty")
write.csv(data.frame(wl = 400:2500, r_forest = res_wire_inf_lib), file.path(outdir, "ref_inform_liberty.csv"), row.names = FALSE)

lut_wire_inf_fcx <- data.frame(c(list(Cab = 40, Car = 10, EWT = 0.015, LMA = 0.01, Cs = 0.1, N = 1.8, fqe = 0.01, Cx = 0.5,
                                       Prot = 0.001, CBC = 0.008, Anth = 1, alpha = 40), inform_common))
res_wire_inf_fcx <- ToolsRTM::inform(inputLUT = lut_wire_inf_fcx, rsoil = rsoil_2101, LeafModel = "Fluspect-B-Cx")
write.csv(data.frame(wl = 400:2400, r_forest = res_wire_inf_fcx), file.path(outdir, "ref_inform_fluspectcx.csv"), row.names = FALSE)

cat("Leaf-model wiring reference export done\n")

cat("SCOPEinR BSM reference export done\n")

## ---- 13. Leaf biochemistry reference (get.biochemical, C3 + C4) ----
## Written directly to the scopeinpython package's own tests/refdata/ (not
## staged in python/scratch/_refdata/ first, unlike the ToolsRTM-side sections
## above) since this data is scopeinpython-only and there's no separate
## "review in _refdata before bundling" step needed for it.
scopeinpython_out_dir <- file.path(root, "python/scopeinpython/tests/refdata")
dir.create(scopeinpython_out_dir, showWarnings = FALSE, recursive = TRUE)

TDP_bio <- SCOPEinR::define_temp_response_biochem(getTDP = TRUE)
data.meteo_bio <- list(Q = 800, Cs = 400, Temp = 25, eb = 15, Oa = 209, p = 1013)
biochem_fields <- c("A","Ci","rcw","gs","RH","Vcmax","Rd","Ja","ps","ps_rel","Kd","Kn","NPQ","Kp",
                     "eta","qE","fs","SIF","fo0","fm0","fo","fm","qQ","Phi_N")

data.leafbio_c3 <- list(Type = "C3", stressfactor = 1, Vcmax25 = 60, BallBerry0 = 0.01, BallBerrySlope = 9,
                         Rdparam = 0.015, Kn0 = 2.48, Knalpha = 2.83, Knbeta = 0.114, TDP = TDP_bio)
data.opts_bio <- data.frame(Value = rep(1, 10), Name = paste0("opt", 1:10))  # row 7 = tempcor = 1
res_bio_c3 <- SCOPEinR::get.biochemical(data.leafbio = data.leafbio_c3, data.meteo = data.meteo_bio,
                                         data.opts = data.opts_bio, fV = 1, get.plots = FALSE)
write.csv(as.data.frame(t(sapply(biochem_fields, function(f) res_bio_c3[[f]]))),
          file.path(scopeinpython_out_dir, "ref_biochemical_c3.csv"), row.names = FALSE)

## NB: Type='C4' with tempcor=0 crashes in R itself ("object 'Vcmax' not
## found") -- biochemical.R only ever assigns Vcmax/Rd for C4 inside the
## tempcor==1 branch, and the separate C3-only `if (Type == "C3")` block is
## the only other place they're set. A real bug in biochemical.R (C4
## requires tempcor=1), reproduced (as the same NameError) rather than
## routed around in the Python port -- see scopeinpython/biochemical.py.
data.leafbio_c4 <- list(Type = "C4", stressfactor = 1, Vcmax25 = 30, BallBerry0 = 0, BallBerrySlope = 4,
                         Rdparam = 0.025, Kn0 = 2.48, Knalpha = 2.83, Knbeta = 0.114, TDP = TDP_bio)
res_bio_c4 <- SCOPEinR::get.biochemical(data.leafbio = data.leafbio_c4, data.meteo = data.meteo_bio,
                                         data.opts = data.opts_bio, fV = 1, get.plots = FALSE)
write.csv(as.data.frame(t(sapply(biochem_fields, function(f) res_bio_c4[[f]]))),
          file.path(scopeinpython_out_dir, "ref_biochemical_c4.csv"), row.names = FALSE)

cat("Leaf biochemistry reference export done\n")

## ---- 14. SCOPE-variant Fluspect-B-Cx reference (getFluspect.Cx.SCOPE,
## the function SCOPE's own leaf-optics pipeline actually calls -- genuinely
## different from ToolsRTM::getFluspect.Cx despite similar code, see
## python/scopeinpython/src/scopeinpython/fluspect.py) ----
op_scope <- SCOPEinR::optipar2021.Pro.CX  # the real default fluspect_mSCOPE.R uses (optipar2017.ProspectD lacks Kp/Kcbc data)
op_scope_cols <- c("wl","nr","Kab","Kca","Ks","Kw","Kdm","KcaV","KcaZ","Kant","Kp","Kcbc","phi","phiI","phiII")
write.csv(as.data.frame(op_scope[op_scope_cols]),
          file.path(root, "python/scopeinpython/src/scopeinpython/data/optipar_fluspect_scope.csv"), row.names = FALSE)

lut_flus_scope <- data.frame(Cab = 40, Car = 10, EWT = 0.015, LMA = 0.01, Cs = 0.1, N = 1.8, fqe = 0.01, Cx = 0.5,
                              Prot = 0.001, CBC = 0.008, Anth = 1)
for (step_i in c(5, 1)) {
  res_fs <- SCOPEinR::getFluspect.Cx.SCOPE(inputsLeaf = lut_flus_scope, inputsOptipar = op_scope,
                                            version = 'Cx', step = step_i)
  label <- paste0("step", step_i)
  write.csv(data.frame(wl = res_fs$lambda, refl = res_fs$refl, tran = res_fs$tran,
                        kChlrel = res_fs$kChlrel, kCarrel = res_fs$kCarrel),
            file.path(scopeinpython_out_dir, paste0("ref_fluspect_scope_", label, "_refltran.csv")), row.names = FALSE)
  write.csv(res_fs$Mb, file.path(scopeinpython_out_dir, paste0("ref_fluspect_scope_", label, "_Mb.csv")), row.names = FALSE)
  write.csv(res_fs$Mf, file.path(scopeinpython_out_dir, paste0("ref_fluspect_scope_", label, "_Mf.csv")), row.names = FALSE)
}

cat("SCOPE-variant Fluspect-Cx reference export done\n")

## ---- 15. fluspect_mSCOPE reference (get.fluspect_mSCOPE, the multi-layer
## (mSCOPE) wrapper around getFluspect.Cx.SCOPE -- computes leaf optics once
## per distinct leaf-biochemistry profile layer, then replicates across the
## canopy sublayers each profile layer spans, see
## python/scopeinpython/src/scopeinpython/fluspect_mscope.py) ----
spectral_full <- SCOPEinR::get.spectra.SCOPE(getSpectral = TRUE)
mly <- list(nly = 3, pLAI = c(0.5, 1.0, 0.5),
            pCab = c(50, 35, 20), pEWT = c(0.02, 0.015, 0.01),
            pCar = c(12, 9, 6), pLMA = c(0.012, 0.01, 0.008),
            pCs = c(0.05, 0.1, 0.15), pN = c(2.0, 1.8, 1.5))
leafbio_ml <- list(Cx = 0.5, fqe = 0.01, Prot = 0.001, CBC = 0.008, Anth = 1)
nl_ml <- 10
res_ml <- SCOPEinR::get.fluspect_mSCOPE(mly = mly, spectral = spectral_full, leafbio = leafbio_ml,
                                         soil = list(rs_thermal = 0.1), optipar = NULL, nl = nl_ml,
                                         step = 5, get.plots = FALSE)
write.csv(res_ml$refl, file.path(scopeinpython_out_dir, "ref_fluspect_mscope_refl.csv"), row.names = FALSE)
write.csv(res_ml$tran, file.path(scopeinpython_out_dir, "ref_fluspect_mscope_tran.csv"), row.names = FALSE)
write.csv(res_ml$kChlrel, file.path(scopeinpython_out_dir, "ref_fluspect_mscope_kChlrel.csv"), row.names = FALSE)
write.csv(res_ml$kCarrel, file.path(scopeinpython_out_dir, "ref_fluspect_mscope_kCarrel.csv"), row.names = FALSE)
write.csv(matrix(res_ml$Mb, nrow = dim(res_ml$Mb)[1] * dim(res_ml$Mb)[2], ncol = dim(res_ml$Mb)[3]),
          file.path(scopeinpython_out_dir, "ref_fluspect_mscope_Mb.csv"), row.names = FALSE)
write.csv(matrix(res_ml$Mf, nrow = dim(res_ml$Mf)[1] * dim(res_ml$Mf)[2], ncol = dim(res_ml$Mf)[3]),
          file.path(scopeinpython_out_dir, "ref_fluspect_mscope_Mf.csv"), row.names = FALSE)
write.csv(data.frame(phiI = res_ml$phiI, phiII = res_ml$phiII),
          file.path(scopeinpython_out_dir, "ref_fluspect_mscope_phi.csv"), row.names = FALSE)
cat("fluspect_mSCOPE dims: Mb=", paste(dim(res_ml$Mb), collapse="x"),
    " refl=", paste(dim(res_ml$refl), collapse="x"), "\n")

cat("fluspect_mSCOPE reference export done\n")
