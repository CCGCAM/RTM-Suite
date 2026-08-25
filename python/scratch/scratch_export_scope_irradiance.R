# Exports SCOPEinR's bundled default irradiance (Esun_[[1]]/Esky_[[1]]),
# used by get.SCOPE() when options.irradiance == 0, as package data for the
# Python port (scopeinpython.data.default_irradiance.csv). Run with
# devtools::load_all() against current R source, not library() against a
# possibly-stale install.

devtools::load_all("SCOPEinR")

wl <- SCOPEinR::get.spectra.SCOPE(getSpectral = TRUE)[["wlIrrad"]]
Esun_ <- SCOPEinR::Esun_[[1]]
Esky_ <- SCOPEinR::Esky_[[1]]

stopifnot(length(wl) == length(Esun_), length(wl) == length(Esky_))

out <- data.frame(wave = wl, Esun_ = Esun_, Esky_ = Esky_)
write.csv(out, "python/scopeinpython/src/scopeinpython/data/default_irradiance.csv", row.names = FALSE)
cat("wrote", nrow(out), "rows\n")
