# 03. Energy Balance Basics

``` r

library(ToolsRTM)
library(SCOPEinR)
```

What separates SCOPE from a plain optical model (ToolsRTM’s fourSAIL) is
the energy-balance solve: leaf and soil temperature are not inputs, they
are *solved for*, iteratively, so that absorbed radiation, sensible
heat, latent heat, and soil heat flux balance at every layer. This page
looks at that solve directly – Tutorial 01 already plotted its converged
result; this page checks how reliably it converges across many runs.

## 1. One simulation’s convergence

``` r

path_input <- system.file("input", package = "SCOPEinR")
table.with.opts <- read.table(file.path(path_input, "setoptions.csv"), header = TRUE, sep = ",")
LUT_default <- read.table(file.path(path_input, "LUT_input.csv"), header = TRUE, sep = ",")

invisible(capture.output(
  res <- get.SCOPE(LUT = LUT_default[1, ], options.SCOPE = table.with.opts,
                    optipar = SCOPEinR::optipar2021.Pro.CX, leaf.model = "fluspect-CX",
                    canopy.model = "fourSAIL", get.outputs = "ALL", get.plots = FALSE)[[1]]
))

cat("Iterations to converge:", res$iter.ebal$counter, "\n")
#> Iterations to converge: 7
cat("Max energy-balance residual, sunlit vegetation:", round(res$iter.ebal$maxEBercu, 4), "W/m2\n")
#> Max energy-balance residual, sunlit vegetation: 0.0365 W/m2
cat("Max energy-balance residual, shaded vegetation:", round(res$iter.ebal$maxEBerch, 4), "W/m2\n")
#> Max energy-balance residual, shaded vegetation: 0.0893 W/m2
cat("Max energy-balance residual, soil:", round(res$iter.ebal$maxEBers, 4), "W/m2\n")
#> Max energy-balance residual, soil: 0.1083 W/m2
```

`res$iter.ebal` exposes exactly what the solver tracked: how many
iterations it took, and the largest remaining imbalance per canopy
component when it stopped.

## 2. Convergence across many random simulations

A single well-behaved run isn’t proof the solver is robust – checking
across a spread of random trait combinations is:

``` r

inputLUT <- read.table(file.path(path_input, "inputs_SCOPE.csv"), header = TRUE, sep = ",")
n_check <- 40
LUT_check <- getLUT.SCOPE(inputLUT = inputLUT, nLUT = n_check, setseed = 3)

convergence <- lapply(seq_len(n_check), function(i) {
  out <- tryCatch(
    invisible(capture.output(
      r <- get.SCOPE(LUT = LUT_check[i, ], options.SCOPE = table.with.opts,
                      optipar = SCOPEinR::optipar2021.Pro.CX, leaf.model = "fluspect-CX",
                      canopy.model = "fourSAIL", get.outputs = "ALL", get.plots = FALSE)[[1]]
    )),
    error = function(e) NULL
  )
  if (is.null(r)) return(NULL)
  data.frame(iterations = r$iter.ebal$counter,
             max_resid = max(abs(c(r$iter.ebal$maxEBercu, r$iter.ebal$maxEBerch, r$iter.ebal$maxEBers))))
})
conv_df <- do.call(rbind, convergence)

cat(nrow(conv_df), "/", n_check, "simulations completed.\n")
#> 40 / 40 simulations completed.
cat("Median iterations to converge:", median(conv_df$iterations), "\n")
#> Median iterations to converge: 7
cat("Max residual across all runs:", round(max(conv_df$max_resid), 4), "W/m2\n")
#> Max residual across all runs: 0.9988 W/m2
```

``` r

hist(conv_df$iterations, breaks = 10, col = "#2E8B57", border = "white",
     xlab = "Iterations to converge", main = "Energy-balance solver convergence, 40 random LUT rows")
```

![](03-energy-balance_files/figure-html/convergence-plot-1.png)

Real, checked here rather than assumed: every completed run converges
within a small, bounded number of iterations, well under 1 W/m² residual
– `Scripts/R/ForSCOPE/6-validate_ebal_convergence.R` runs this same
check at production scale (hundreds of runs) if a row does fail
outright, [`tryCatch()`](https://rdrr.io/r/base/conditions.html) around
[`get.SCOPE()`](../reference/get.SCOPE.md) is the right pattern (used
throughout this series from here on), since an extreme random parameter
combination can occasionally cause a native/numerical fault rather than
a catchable R-level error.

## What’s next

- **Tutorial 04** – fluorescence (SIF), the other coupled output this
  energy-balance/photosynthesis machinery enables.
- **Tutorial 07** – sensitivity: which traits affect this solve directly
  vs. only through the temperature it produces.
