#' latin_hypercube_input function
#'
#' \code{latin_hypercube_input} draws a maximin Latin Hypercube sample of
#' \code{n_spectra} parameter sets over the parameter ranges given in
#' \code{tab}, and writes them to a CSV file (\code{lh_ts.csv}) for use as
#' SCOPE simulation input (e.g. to generate a synthetic training/verification
#' dataset). If \code{LIDFa} and \code{LIDFb} are among the sampled
#' variables, the sampled pair is transformed so that \code{abs(LIDFa +
#' LIDFb) <= 1} is respected.
#'
#' @param tab data.frame. Parameter bounds table with columns \code{include} (logical, whether the variable is sampled), \code{lower}/\code{upper} (parameter bounds), and \code{variable} (parameter name). Defaults to reading \code{input/dataset for_verification/input_borders.csv}.
#' @param n_spectra integer. Number of parameter sets (LHS samples) to draw. Default 30.
#' @param outdir character. Output directory in which \code{lh_ts.csv} is written (and, if \code{tab} is missing, from which \code{input_borders.csv} is read). Default \code{input/dataset for_verification}.
#'
#' @return No return value; called for its side effect of writing the sampled parameter table to \code{file.path(outdir, 'lh_ts.csv')} and printing a summary message. Errors if that file already exists.
#' @export
#'
#' @author 	Christiaan van der Tol(Original version in Matlab)
#' @author Carlos Camino (Ported version into R)
#'
#' @examples
#' \dontrun{
#' latin_hypercube_input(tab, n_spectra = 30, outdir = "input/dataset for_verification")
#' }
latin_hypercube_input <- function(tab = read.table(file.path('input', 'dataset for_verification', 'input_borders.csv'), header = TRUE),
                                  n_spectra = 30, outdir = file.path('input', 'dataset for_verification')) {

  if (missing(tab)) {
    tab <- read.table(file.path(outdir, 'input_borders.csv'), header = TRUE)
    n_spectra <- 30
  }

  out_file <- file.path(outdir, 'lh_ts.csv')
  if (file.exists(out_file)) stop(sprintf('`%s` file already exists, delete it first', out_file))

  include <- as.logical(tab$include)
  lb <- as.numeric(tab$lower[include])
  ub <- as.numeric(tab$upper[include])
  varnames <- as.character(tab$variable[include])

  # one row - one set of parameters
  lh <- lhs::maximinLHS(n = n_spectra, k = sum(include))
  params <- t(t((ub-lb) * lh) + lb)

  if ('LIDFa' %in% varnames) {
    # abs(LIDFa + LIDFb) <= 1
    i_lidfa <- which(varnames == 'LIDFa')
    i_lidfb <- which(varnames == 'LIDFb')
    lidfa <- params[, i_lidfa]
    lidfb <- params[, i_lidfb]
    params[, i_lidfa] <- (lidfa + lidfb) / 2
    params[, i_lidfb] <- (lidfa - lidfb) / 2
  }

  t <- data.frame(params)
  colnames(t) <- varnames
  t$t <- seq_len(nrow(t))
  write.table(t, out_file, row.names = FALSE)

  varnames_in <- paste(varnames, collapse = ', ')
  cat(sprintf('Sampled %i parameters: %s\n', length(varnames), varnames_in))
  cat(sprintf('Saved lut input (parameters) in `%s`\n', out_file))
}
