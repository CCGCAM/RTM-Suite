## One-off: verify get.marmit.rsoil()'s db_root argument loads all 8 official
## MARMIT databases from the repo's own databases/ folder correctly (id=1 in
## each, real files -- not synthetic). Run from the repo root:
##   Rscript python/scratch/scratch_verify_marmit_databases.R
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
root <- normalizePath(file.path(dirname(this_file), "..", ".."))

devtools::load_all(file.path(root, "ToolsRTM/R"), quiet = TRUE)

db_root <- file.path(root, "databases")
databases <- c("Bablet_2016", "Dupiau_2020", "Humper_2015", "Lesaignoux_2008",
               "Liu_2002", "Lobell_2002", "Marcq_2012", "Philpot_2014")

cat("=== db_root =", db_root, "===\n")
for (db in databases) {
  soil <- get.marmit.rsoil(database = db, id = 1, L = 0.05, eps = 0.3, db_root = db_root)
  stopifnot(all(is.finite(soil$rsoil.wet)), all(is.finite(soil$rsoil.dry)))
  stopifnot(all(soil$rsoil.wet <= soil$rsoil.dry + 1e-9))
  stopifnot(soil$SMC >= 0, soil$SMC <= 100)
  cat(sprintf("%-16s OK  SMC = %.2f%%\n", db, soil$SMC))
}

## Bundled Bablet_2016 (no db_root) must match the same database via db_root
bundled <- get.marmit.rsoil(database = "Bablet_2016", id = 1, L = 0.05, eps = 0.3)
external <- get.marmit.rsoil(database = "Bablet_2016", id = 1, L = 0.05, eps = 0.3, db_root = db_root)
stopifnot(isTRUE(all.equal(bundled$SMC, external$SMC, tolerance = 1e-6)))
stopifnot(isTRUE(all.equal(bundled$rsoil.wet, external$rsoil.wet, tolerance = 1e-6)))
cat("Bundled Bablet_2016 == databases/Bablet_2016 via db_root: OK\n")

## Unknown database error path
err <- tryCatch(get.marmit.rsoil(database = "Nonexistent", db_root = db_root),
                error = function(e) conditionMessage(e))
stopifnot(grepl("Bablet_2016", err))
cat("Unknown-database error message lists available databases: OK\n")

cat("\nALL 8 MARMIT DATABASES VERIFIED\n")
