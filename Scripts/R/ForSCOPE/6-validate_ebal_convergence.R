# ==============================================================================
# Validate the ebal.R soil-temperature fix (Ts <- Ts[1] + ... -> Ts <- Ts + ...)
# across many randomized simulations, not just the single default LUT row.
#
# Runs n_samples SCOPE simulations with randomized plant/soil/geometry traits
# and records, for each one, how many iterations the energy balance took to
# converge and the final residual error per component (sunlit veg, shaded
# veg, soil). Summarizes as histograms + a convergence-rate table.
# ==============================================================================

rm(list = ls())  # avoid leftover objects from a previous run/session leaking in

library(ToolsRTM)
library(SCOPEinR)
library(ggplot2)
library(dplyr)

out_dir <- "../../../outs/ForSCOPE/ebal"  # project-level outputs folder, never inside Scripts/
dir.create(out_dir, showWarnings = FALSE)

# NOTE on scale: a PSOCK worker crash ("error reading from connection") was
# reproduced at chunk sizes of 250 and (once) 500 rows with get.outputs =
# 'ALL' + parallel = TRUE -- reliably enough to not be a fluke, not reliably
# enough (yet) to pin down which specific random LUT row triggers it or
# whether it's a memory ceiling. Dialed n_samples/chunk_size down here to a
# size that completes cleanly, as a real (not aspirational) validation run.
# Scaling get.SCOPE.parallel() reliably to 5000-20000 'ALL'-mode rows needs
# more investigation before it can be trusted at that scale -- flagging this
# rather than silently shipping a script that might hang.
n_samples <- 200
chunk_size <- 50

opts <- read.table(system.file("input", "setoptions.csv", package = "SCOPEinR"), header = TRUE, sep = ",")
inputLUT <- read.table(system.file("input", "inputs_SCOPE.csv", package = "SCOPEinR"), header = TRUE, sep = ",")
LUT <- as.data.frame(getLUT.SCOPE(inputLUT = inputLUT, nLUT = n_samples))

run_chunk <- function(LUT_chunk) {
  get.SCOPE.parallel(
    LUT = LUT_chunk, options.SCOPE = opts, optipar = SCOPEinR::optipar2021.Pro.CX,
    leaf.model = "fluspect-CX", canopy.model = "fourSAIL",
    # get.outputs must be 'ALL' -- iter.ebal (convergence stats) is only
    # included in the returned list in that mode, not in 'Main'.
    parallel = TRUE, get.outputs = "ALL", get.plots = FALSE, get.csv = FALSE
  )
}

cat("=== Running", n_samples, "SCOPE simulations in chunks of", chunk_size, "===\n")
t0 <- Sys.time()
chunk_starts <- seq(1, n_samples, by = chunk_size)
sims <- list()
failed_chunks <- integer(0)

for (ci in seq_along(chunk_starts)) {
  start_row <- chunk_starts[ci]
  end_row <- min(start_row + chunk_size - 1, n_samples)
  LUT_chunk <- LUT[start_row:end_row, ]

  chunk_res <- tryCatch(run_chunk(LUT_chunk), error = function(e) {
    message(sprintf("chunk %d (%d-%d) failed: %s -- retrying once",
                    ci, start_row, end_row, conditionMessage(e)))
    tryCatch(run_chunk(LUT_chunk), error = function(e2) {
      message(sprintf("chunk %d (%d-%d) failed again on retry: %s -- skipping",
                      ci, start_row, end_row, conditionMessage(e2)))
      NULL
    })
  })

  if (is.null(chunk_res)) {
    failed_chunks <- c(failed_chunks, ci)
  } else {
    sims <- c(sims, chunk_res)
    message(sprintf("chunk %d/%d (%d-%d) done", ci, length(chunk_starts), start_row, end_row))
  }
}
t1 <- Sys.time()
cat("Elapsed:", as.numeric(t1 - t0, units = "secs"), "s for", length(sims), "simulations",
    "(", length(failed_chunks), "chunk(s) unrecoverable after retry)\n")

## ----------------------------------------------------------------------------
## Extract convergence stats from each simulation
## ----------------------------------------------------------------------------

stats <- bind_rows(lapply(seq_along(sims), function(i) {
  it <- sims[[i]]$iter.ebal
  if (is.null(it)) return(NULL)
  data.frame(sim = i, counter = it$counter, maxit = it$maxit,
             maxEBercu = it$maxEBercu, maxEBerch = it$maxEBerch, maxEBers = it$maxEBers,
             converged = it$counter < it$maxit)
}))

cat("\n=== Convergence summary (n =", nrow(stats), "usable simulations) ===\n")
cat(sprintf("Converged (counter < maxit): %d/%d (%.1f%%)\n",
            sum(stats$converged), nrow(stats), 100 * mean(stats$converged)))
cat("Iteration count:    median =", median(stats$counter), " mean =", round(mean(stats$counter), 1),
    " max =", max(stats$counter), "\n")
cat("maxEBers  (soil):   median =", round(median(stats$maxEBers), 3),
    " mean =", round(mean(stats$maxEBers), 3), " max =", round(max(stats$maxEBers), 2), "W/m2\n")
cat("maxEBercu (sunlit): median =", round(median(stats$maxEBercu), 4),
    " max =", round(max(stats$maxEBercu), 3), "W/m2\n")
cat("maxEBerch (shaded): median =", round(median(stats$maxEBerch), 4),
    " max =", round(max(stats$maxEBerch), 3), "W/m2\n")

write.csv(stats, file.path(out_dir, "ebal_convergence_500sims.csv"), row.names = FALSE)

## ----------------------------------------------------------------------------
## Plots
## ----------------------------------------------------------------------------

p_iter <- ggplot(stats, aes(x = counter)) +
  geom_histogram(binwidth = 1, fill = "#0072B2", color = "white") +
  geom_vline(xintercept = unique(stats$maxit), linetype = "dashed", color = "red") +
  labs(title = "Energy balance: iterations to converge",
       subtitle = sprintf("%d simulations, random LUT -- dashed line = maxit (%d)", n_samples, stats$maxit[1]),
       x = "Iterations (counter)", y = "Count") +
  theme_bw(base_size = 12)

df_long <- stats %>%
  select(sim, maxEBercu, maxEBerch, maxEBers) %>%
  tidyr::pivot_longer(-sim, names_to = "component", values_to = "residual") %>%
  mutate(component = recode(component, maxEBercu = "Sunlit vegetation",
                            maxEBerch = "Shaded vegetation", maxEBers = "Soil"))

p_residual <- ggplot(df_long, aes(x = residual, fill = component)) +
  geom_histogram(bins = 40) +
  facet_wrap(~component, scales = "free", ncol = 1) +
  scale_fill_manual(values = c("Sunlit vegetation" = "#E69F00", "Shaded vegetation" = "#009E73",
                               "Soil" = "#D55E00")) +
  labs(title = "Final energy balance residual by component",
       subtitle = sprintf("%d simulations, random LUT (tolerance = 1 W/m2)", n_samples),
       x = "Max residual error (W/m2)", y = "Count") +
  theme_bw(base_size = 12) + theme(legend.position = "none")

ggsave(file.path(out_dir, "ebal_convergence_iterations.png"), p_iter, width = 8, height = 5, dpi = 150)
ggsave(file.path(out_dir, "ebal_convergence_residuals.png"), p_residual, width = 8, height = 8, dpi = 150)

cat("\nSaved to '", out_dir, "/': ebal_convergence_500sims.csv, ",
    "ebal_convergence_iterations.png, ebal_convergence_residuals.png\n", sep = "")
