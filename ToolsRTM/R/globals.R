# Suppress R CMD check NOTEs ("no visible binding for global variable")
# for non-standard-evaluation data-frame/ggplot columns, foreach/apply loop
# variables, and the package's own lazy-loaded datasets referenced without
# an explicit `ToolsRTM::` prefix from inside the package's own functions.
utils::globalVariables(c(
  "Eo", "Extraterrestrial_irradiance", "I.Johnson", "ID", "IncludeModel",
  "Kant", "Kca", "Kcab", "Kcbc", "Kp", "LUT", "RFL", "Response", "Sensor",
  "Trait", "average", "band", "data", "dataSpec_PDB", "dataSpec_PRO",
  "group", "i", "ind", "inputs", "inputs.SPART", "irrad", "leafbio", "lm",
  "median", "n_bands.points", "na.pass", "percentile_10", "percentile_25",
  "percentile_75", "percentile_90", "pred", "qchisq", "qgamma", "qnbinom",
  "qpois", "quantile", "quantile.25", "quantile.75", "qweibull", "r_leaf",
  "r_understorey", "rad.toa", "rast.2020", "raster", "reflectance",
  "rfl.soil", "rfl.toa", "rfl.toc", "rfl.toc.brdf", "se", "spectra_band",
  "stack.image", "t_leaf", "value", "values", "variable", "wave",
  "wavelength", "wl", "y"
))
