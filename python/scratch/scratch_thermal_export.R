## Reference export for the thermal.py building blocks (Monin-Obukhov,
## resistances, heatfluxes) -- the first pieces of the SCOPE thermal
## energy-balance chain (task 5).
## Repo root, derived from this script's own location (python/scratch/<file>.R)
## rather than hardcoded, so it works regardless of who runs it or from where.
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
root <- normalizePath(file.path(dirname(this_file), "..", ".."))
outdir <- file.path(root, "python/scopeinpython/tests/refdata")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

suppressMessages(devtools::load_all(file.path(root, "SCOPEinR/R"), quiet = TRUE))

## ---- Monin-Obukhov, a few representative (ustar, Ta, H) triples ----
ustar_v <- c(0.3, 0.15, 0.6)
Ta_v <- c(20, 25, 15)
H_v <- c(50, -20, 150)
L_v <- SCOPEinR::get.Monin.Obukhov(data.meteo = list(ustar = ustar_v, Ta = Ta_v), H = H_v)
write.csv(data.frame(ustar = ustar_v, Ta = Ta_v, H = H_v, L = L_v),
          file.path(outdir, "ref_monin_obukhov.csv"), row.names = FALSE)

## ---- resistances, one representative case each for unstable/stable/neutral L ----
res_cases <- list(
  unstable = list(L = -50, z = 40, u = 3, d = 12, hc = 15, z0m = 1.5, LAI = 3, Cd = 0.3, rwc = 100, leafwidth = 0.1, rbs = 50),
  stable   = list(L = 50,  z = 40, u = 2, d = 12, hc = 15, z0m = 1.5, LAI = 3, Cd = 0.3, rwc = 100, leafwidth = 0.1, rbs = 50),
  neutral  = list(L = -1e6, z = 40, u = 3, d = 12, hc = 15, z0m = 1.5, LAI = 3, Cd = 0.3, rwc = 100, leafwidth = 0.1, rbs = 50)
)
res_rows <- list()
for (nm in names(res_cases)) {
  cs <- res_cases[[nm]]
  data.soil <- list(rbs = cs$rbs)
  data.canopy <- list(Cd = cs$Cd, LAI = cs$LAI, rwc = cs$rwc, zo = cs$z0m, d = cs$d, hc = cs$hc, leafwidth = cs$leafwidth)
  data.meteo <- list(z = cs$z, u = cs$u, L = cs$L)
  out <- SCOPEinR::get.resistances(data.soil, data.canopy, data.meteo)
  res_rows[[nm]] <- data.frame(case = nm, L = cs$L, ustar = out$ustar, uz0 = out$uz0, Kh = out$Kh,
                                rai = out$rai, rar = out$rar, rac = out$rac, rws = out$rws,
                                raa = out$raa, rawc = out$rawc, raws = out$raws)
}
write.csv(do.call(rbind, res_rows), file.path(outdir, "ref_resistances.csv"), row.names = FALSE)

## ---- heatfluxes, vectorised over a few (Tc,ea,Ta,Ca,Ci) combinations ----
ra_v <- c(50, 60, 45)
rs_v <- c(200, 250, 180)
Tc_v <- c(22, 28, 18)
ea_v <- c(15, 18, 12)
Ta_v2 <- c(20, 25, 16)
Ca_v <- c(400, 400, 400)
Ci_v <- c(280, 300, 260)
e_to_q <- 0.622 / 1013
hf <- SCOPEinR::get.heatfluxes(ra_v, rs_v, Tc_v, ea_v, Ta_v2, e_to_q, Ca_v, Ci_v)
write.csv(data.frame(ra = ra_v, rs = rs_v, Tc = Tc_v, ea = ea_v, Ta = Ta_v2, Ca = Ca_v, Ci = Ci_v,
                      lambda_ = hf$lambda, s = hf$s, lE = hf$lE, H = hf$H, ec = hf$ec, Cc = hf$Cc),
          file.path(outdir, "ref_heatfluxes.csv"), row.names = FALSE)

cat("Thermal building-blocks reference export done.\n")
