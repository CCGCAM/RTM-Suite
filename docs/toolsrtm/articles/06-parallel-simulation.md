# 06. Large-Scale and Parallel RTM Simulation

``` r

library(ToolsRTM)
library(doParallel)
library(foreach)
```

Every simulation so far used a plain
[`sapply()`](https://rdrr.io/r/base/lapply.html)/loop over LUT rows –
fine for hundreds of rows, slow for thousands. This page shows the
`doParallel`/`foreach` pattern this package’s own course scripts use
(`Scripts/R/ForPROSAIL/1-GetSimulationsLUTs.R`) to distribute
simulations across CPU cores.

## 1. Build a LUT

``` r

n.samples <- 300
LUT <- as.data.frame(getLUT(inputs = ToolsRTM::inputsPROSAIL, nLUT = n.samples, setseed = 1234))
rsoil <- rep(0.15, 2101)
```

## 2. Sequential baseline

``` r

t_seq <- system.time({
  sims_seq <- lapply(seq_len(n.samples), function(i) {
    sail_i <- foursail(inputLUT = LUT[i, ], rsoil = rsoil, LeafModel = "PROSPECT-D")
    sail_i$rsot
  })
})
cat("Sequential:", round(t_seq[["elapsed"]], 2), "s for", n.samples, "simulations\n")
#> Sequential: 1.55 s for 300 simulations
```

## 3. Parallel with `doParallel`/`foreach`

``` r

no_cores <- max(1, parallel::detectCores() - 2)
cl <- makeCluster(no_cores)
registerDoParallel(cl)

t_par <- system.time({
  sims_par <- foreach(i = 1:n.samples, .packages = "ToolsRTM") %dopar% {
    sail_i <- foursail(inputLUT = LUT[i, ], rsoil = rsoil, LeafModel = "PROSPECT-D")
    sail_i$rsot
  }
})
stopCluster(cl)

cat("Parallel (", no_cores, "cores):", round(t_par[["elapsed"]], 2), "s for", n.samples, "simulations\n")
#> Parallel ( 10 cores): 3.11 s for 300 simulations
cat("Speedup:", round(t_seq[["elapsed"]] / t_par[["elapsed"]], 2), "x\n")
#> Speedup: 0.5 x
```

Three things matter for this pattern to actually work, not just look
parallel:

- **[`makeCluster()`](https://rdrr.io/r/parallel/makeCluster.html)/[`registerDoParallel()`](https://rdrr.io/pkg/doParallel/man/registerDoParallel.html)/[`stopCluster()`](https://rdrr.io/r/parallel/makeCluster.html)**
  always paired – an unclosed cluster leaks worker processes across the
  whole R session, compounding if this chunk runs more than once
  (e.g. re-knitting interactively).
- **`.packages = "ToolsRTM"`** inside
  [`foreach()`](https://rdrr.io/pkg/foreach/man/foreach.html) – each
  worker process starts fresh, with no packages attached; without this,
  [`foursail()`](../reference/foursail.md) isn’t found on the workers
  and every iteration errors.
- **`no_cores <- detectCores() - 2`**, not
  [`detectCores()`](https://rdrr.io/r/parallel/detectCores.html) itself
  – leaves headroom for the OS and this R session’s own main process, a
  convention this package’s own course scripts (`ForPROSAIL`,
  `Sensibility`) use throughout.

``` r

identical(sims_seq[[1]], sims_par[[1]])  # same result, sequential or parallel
#> [1] TRUE
```

## 4. When parallel is (and isn’t) worth it

``` r

cat("Per-simulation cost: sequential", round(1000 * t_seq[["elapsed"]] / n.samples, 2),
    "ms, parallel", round(1000 * t_par[["elapsed"]] / n.samples, 2), "ms\n")
#> Per-simulation cost: sequential 5.17 ms, parallel 10.37 ms
```

[`foursail()`](../reference/foursail.md) itself is fast (milliseconds
per call) – the parallel speedup above has to overcome the fixed
overhead of starting a cluster and shipping `ToolsRTM`/data to every
worker, so it only pays off once `n.samples` is large enough (thousands,
not hundreds) or each simulation is individually expensive
(e.g. [`SPART()`](../reference/SPART.md)’s atmosphere calculations,
Tutorial 03; the full SCOPE energy balance, covered in the SCOPEinR
package’s own tutorials). For a few hundred rows, the sequential
[`lapply()`](https://rdrr.io/r/base/lapply.html)/[`sapply()`](https://rdrr.io/r/base/lapply.html)
used throughout Tutorials 01-05 is often simpler and not meaningfully
slower once cluster startup is accounted for.

## 5. Filtering non-finite rows

Large parallel runs sometimes include a few unstable simulations (e.g.
Liberty’s leaf solver, Tutorial 02, near the edges of its own parameter
range) that come back as `NaN` rather than erroring outright. Always
check before using the results downstream:

``` r

spectra_mat <- do.call(rbind, sims_par)
ok <- apply(spectra_mat, 1, function(r) all(is.finite(r)))
cat(sum(!ok), "of", n.samples, "simulations were non-finite and would be dropped.\n")
#> 0 of 300 simulations were non-finite and would be dropped.
spectra_mat <- spectra_mat[ok, , drop = FALSE]
LUT_ok <- LUT[ok, ]
```

## What’s next

- **Tutorial 07** – convolving these simulated spectra onto real sensor
  bands.
- **Tutorial 13** – the same simulate-in-parallel pattern as one step in
  a full end-to-end pipeline.
