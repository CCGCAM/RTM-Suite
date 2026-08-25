## One-off: run ToolsRTM's carspls() on the same synthetic dataset used by
## python/toolsrtm/tests/test_inversion.py and save its output for a
## Python-vs-R regression test. Run from the repo root via a local Rscript install:
##   Rscript python/scratch/scratch_export_carspls.R
devtools::load_all("ToolsRTM", quiet = TRUE)

df <- read.csv("python/toolsrtm/tests/refdata/carspls_input.csv")
X <- as.matrix(df[, paste0("V", 1:15)])
y <- df$y

set.seed(1)
res <- carspls(X, y, nLV = 3, fold = 5, scale.pretreat = 1, iteration = 15,
               PartitionType = "interleaved")

write.csv(data.frame(res$Coef), "python/toolsrtm/tests/refdata/carspls_coef.csv", row.names = FALSE)
write.csv(data.frame(RMSECV = res$RMSECV, NumLV = res$NumLV, Nvar = res$Nvar),
          "python/toolsrtm/tests/refdata/carspls_path.csv", row.names = FALSE)
write.csv(data.frame(SelectedVariables = res$SelectedVariables,
                      OptimalIteration = res$Optimal.iteration,
                      MinError = res$MinError),
          "python/toolsrtm/tests/refdata/carspls_selected.csv", row.names = FALSE)
cat("Selected variables:", res$SelectedVariables, "\n")
cat("Optimal iteration:", res$Optimal.iteration, "\n")
cat("Min error:", res$MinError, "\n")
