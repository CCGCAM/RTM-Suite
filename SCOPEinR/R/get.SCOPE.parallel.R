#' Run SCOPE simulations in parallell
#'
#' @param LUT Leaf, biochemistry, viewing angles, meteo, and canopy properties needed for running SCOPE model.
#' @param options.SCOPE Optical leaf properties, the total irradiance if a specific value is provided instead of the usual Modtran output.
#' @param optipar list. Leaf-level optical parameters passed on to \code{\link{get.SCOPE}} for each LUT chunk (e.g. \code{SCOPEinR::optipar2021.Pro.CX}).
#' @param path.out Folder for saving the SCOPE outputs.
#' @param parallel Logical, indicating whether to use parallel processing. Default is TRUE.
#' @param canopy.model Selection of canopy model options available are 'fourSAIL', 'fourSAIL' and 'INFORM'. By default, the fourSAIL model will be used.
#' @param leaf.model Selection of canopy model options available are 'fluspect-CX', 'fluspect-B', 'PROSPECT', 'Liberty'. By default, the fluspect-CX model will be used.
#' @param get.outputs If get.outputs = 'ALL', all variables will be retrieved; if get.outputs = 'Main', only the main variables will be retrieved.
#' By default, SCOPE uses get.outputs = 'ALL', for processing a huge LUT, it is recommended to use get.outputs = 'Main'.
#' @param get.plots Logical, indicating whether to plot the intermediate plots. Default is TRUE.
#' @param get.csv logical, indicating whether to save outputs. Default is TRUE.
#' @param n.cores Integer, indicating the number of cores to use. if this parameter is null or missing Detect cores - 2 will be used
#' @return A list of SCOPE simulations.
#'
#' @examples
#' \dontrun{
#' get.SCOPE.parallel(LUT, options.SCOPE = NULL, optipar = NULL,
#'                    path.out = "output", parallel = TRUE, leaf.model = 'fluspect-CX',
#'                    canopy.model = 'fourSAIL', get.outputs = 'Main', get.plots = TRUE)
#' }
#'
#' @export
get.SCOPE.parallel <- function(LUT, options.SCOPE = NULL,
                               optipar = NULL, path.out = NULL, parallel = TRUE,
                               canopy.model = 'fourSAIL', leaf.model = 'fluspect-CX',
                               get.outputs = 'Main', get.plots = TRUE, get.csv= T, n.cores=8) {
  message("Executing SCOPE 2.1. version ...")


  n.LUT = nrow(LUT)

  if (!(get.outputs %in% c("Main", "ALL"))) {
    message("Warning: 'get.outputs' parameter should be either 'Main' or 'ALL'.")
    stop("Execution halted.")
  }
  if (n.LUT > 150) {
    get_plots <- FALSE
    message("Number of rows in LUT exceeds 150. Disabling plot generation.")
  } else{
    get_plots <- get.plots
  }
  # Define function to process each chunk of LUT data
  process_chunk <- function(chunk) {

    # suppressMessages/Warnings: get.SCOPE() prints a startup message plus an
    # energy-balance convergence report per call, which floods the console
    # once this runs across hundreds/thousands of chunks.
    suppressMessages(suppressWarnings(
      SCOPEinR::get.SCOPE(chunk, n.LUT = nrow(chunk), options.SCOPE = options.SCOPE,
                          optipar = optipar, leaf.model = leaf.model,
                          canopy.model = canopy.model, get.outputs = get.outputs,
                          get.plots = F)
    ))
  }

  if (parallel || is.null(parallel))  {
    # Split LUT into chunks for parallel processing
    # Determine number of CPU cores for parallel processing

    if (missing(n.cores) | is.null(n.cores)){
      num_cores = parallel::detectCores() - 1
    } else {
      num_cores = n.cores
    }
    #num_cores <- ifelse(is.null(server) | server == T, parallel::detectCores() / 10, parallel::detectCores() - 1)
    #rows_to_select <- sample(nrow(LUT), n.LUT) this is not good
    rows_to_select <- c(1:n.LUT)
    chunks <- split(LUT[rows_to_select, ], sort(rep_len(1:num_cores, length.out = n.LUT)))

    # Initialize progress bar
    pb <- progress::progress_bar$new(total = length(chunks), format = "[:bar] :percent eta: :eta")

    # Initialize parallel processing
    cl <- parallel::makeCluster(num_cores)
    doParallel::registerDoParallel(cl)

    # Parallel execution
    start_time <- Sys.time()

    sims <- foreach::foreach(chunk = chunks, .combine = c) %dopar% {
      pb$tick()  # Increment progress bar
      process_chunk(chunk)
    }
    end_time <- Sys.time()

    # Clean up parallel resources
    parallel::stopCluster(cl)

  } else {
    # Non-parallel execution
    start_time <- Sys.time()
    sims <- process_chunk(LUT)
    end_time <- Sys.time()
  }

  # Print execution time
  #cat("Execution time:", end_time - start_time, "\n")
  # Print total number of simulations and total execution time
  total_simulations <- length(sims)
  total_time <- end_time - start_time
  cat("Total simulations:", total_simulations, "\n")
  cat("Total execution time:", total_time, "\n")

  # Save SCOPE outputs as CSV files if required
  if (get.csv) {
    # Indicate that files will be saved
    message("Saving SCOPE outputs...")

    if (is.null(path.out)) {
      path_outs <- 'outs/'
      message('folder "outs" will be used for the SCOPE outputs')

      if (!dir.exists(path_outs)) {
        dir.create(path_outs)
      }
    } else {
      path_outs <- path.out
      if (!dir.exists(path_outs)) {
        dir.create(path_outs)
      }
    }

    SCOPEinR::get.SCOPE.outputs(data.sim = sims, N.sims = length(sims), LUT = LUT,get.outputs=get.outputs,
                                path.out = path_outs, get.plots = get.plots)
  }

  return(sims)

}

