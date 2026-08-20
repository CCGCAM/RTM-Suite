## Repo root, derived from this script's own location (python/scratch/<file>.R)
## rather than hardcoded, so it works regardless of who runs it or from where.
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
root <- normalizePath(file.path(dirname(this_file), "..", ".."))
t0 <- Sys.time()
shiny::runApp(file.path(root, "AEO-Course/Apps/RTMs"), port = 8931, launch.browser = FALSE, host = "127.0.0.1")
