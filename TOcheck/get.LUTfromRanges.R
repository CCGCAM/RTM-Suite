
#' Get a LUT based on a table with Min and Max ranges
#'
#' @param LUT a data.frame with three columns (input,min, max)
#' @param nLUT Number of LUT samples to generate
#' @param set.seed Random seed for reproducibility
#' @param leaf.model Leaf model to use (e.g., 'PROSPECT-PRO', 'PROSPECT-D', 'Liberty', 'FLUSPECT-Cx')
#' @param canopy.model Canopy model to use (e.g., 'fourSAILH', 'INFORM')
#' @param distribution Distribution type for LUT generation ('uniform' or 'gauss')
#'
#' @return List containing LUT for given inputs
#' @export
#'
#' @examples
#' # Generate LUT with PROSPECT-PRO and fourSAILH models using a Gaussian distribution
#' LUT_example <- get.LUTfromRanges(LUT=LUT.range,nLUT = 1000, setseed = 42,
#'                                 leaf.model = 'PROSPECT-PRO',
#'                                 canopy.model = 'fourSAILH',
#'                                 distribution = 'gauss')
#'
#' # Generate LUT with PROSPECT-D and INFORM models using a Uniform distribution
#' LUT_example_uniform <- get.LUTfromRanges(LUT=LUT.range,nLUT = 500, setseed = 123,
#'                                         leaf.model = 'PROSPECT-D',
#'                                         canopy.model = 'INFORM',
#'                                         distribution = 'uniform')
#'
#'


get.LUTfromRanges <- function(LUT=NULL,nLUT=NULL, setseed=1234,
                              leaf.model ='PROSPECT-PRO', canopy.model='fourSAILH',
                              distribution= 'gauss'){

  # Ensure that LUT is provided
  if (is.null(LUT)) stop("The LUT data must be provided.")
  if (is.null(nLUT)){
    message("nLUT is not defined, by default we use a n= 100.")

    nLUT == 100
  }
  if (is.null(setseed)){
    setseed == 1234
    message("set.seed is not defined, by default we use a n= 100.")
  }
  if (is.null(distribution)){
    setseed == 1234
    message("set.seed is not defined, by default we use a n= 100.")
  }

  if (is.null(distribution)) {
    distribution <- 'Uniform'
    message("No distribution specified. Defaulting to 'Uniform'.")
  }

  LUT_ranges <- LUT

  # Check if the LUT has more than three columns and provide a message
  if (ncol(LUT_ranges) > 3) {
    message("Warning: The LUT has more than the expected 3 columns. Only the first 3 columns ('inputs', 'Min', 'Max') will be used.")
  }

  # Assign proper column names (first three columns)
  LUT_ranges <- LUT_ranges[,1:3]
  colnames(LUT_ranges)<- c('inputs', 'Min', 'Max')

  if ((leaf.model == 'PROSPECT-PRO') & (canopy.model == 'fourSAILH')) {
    inputs.leaf = c("N", "Cab", "Car", "Anth", "Cbrown", "EWT", "LMA", "Prot", "CBC", "alpha")
    inputs.canopy = c("LAI", "TypeLidf", "LIDFa", "LIDFb", "hspot", "tts", "tto", "psi", "psoil")

  } else if ((leaf.model == 'PROSPECT-D') & (canopy.model == 'fourSAILH')) {

    inputs.leaf = c("N", "Cab", "Car", "Anth", "Cbrown", "EWT", "LMA", "alpha")
    inputs.canopy = c("LAI", "TypeLidf", "LIDFa", "LIDFb", "hspot", "tts", "tto", "psi", "psoil")

  } else if ((leaf.model == 'PROSPECT-PRO') & (canopy.model == 'INFORM')){

    inputs.leaf = c("N", "Cab", "Car", "Anth", "Cbrown", "EWT", "LMA", "Prot", "CBC", "alpha")
    inputs.canopy = c("LAI", "TypeLidf", "LIDFa", "LIDFb", "hspot", "tts", "tto", "psi", "psoil",
                      "LAIu", "cd", "sd", "h", "skyl")

  } else if ((leaf.model == 'PROSPECT-D') & (canopy.model == 'INFORM')){

    inputs.leaf = c("N", "Cab", "Car", "Anth", "Cbrown", "EWT", "LMA", "alpha")
    inputs.canopy = c("LAI", "TypeLidf", "LIDFa", "LIDFb", "hspot", "tts", "tto", "psi", "psoil",
                      "LAIu", "cd", "sd", "h", "skyl")

  } else if ((leaf.model == 'Liberty') & (canopy.model == 'fourSAILH')){

    inputs.leaf = c("cell.d", "inter.c", "baseline.abs", "leaf.thick", "albino.abs", "Cab", "EWT", "lign.cell", "Nitrogen")
    inputs.canopy = c("LAI", "TypeLidf", "LIDFa", "LIDFb", "hspot", "tts", "tto", "psi", "psoil")

  } else if ((leaf.model == 'Liberty') & (canopy.model == 'INFORM')){

    inputs.leaf = c("cell.d", "inter.c", "baseline.abs", "leaf.thick", "albino.abs", "Cab", "EWT", "lign.cell", "Nitrogen")
    inputs.canopy = c("LAI", "TypeLidf", "LIDFa", "LIDFb", "hspot", "tts", "tto", "psi", "psoil",
                      "LAIu", "cd", "sd", "h", "skyl")

  }  else if ((leaf.model == 'FLUSPECT-B') & (canopy.model == 'fourSAILH')){


    inputs.leaf = c('N',"fqe", "Cab",  "Car", "Cs", "EWT", "LMA", "Cx","alpha")
    inputs.canopy = c("LAI", "TypeLidf", "LIDFa", "LIDFb", "hspot", "tts", "tto", "psi", "psoil")


  } else if ((leaf.model == 'FLUSPECT-B') & (canopy.model == 'INFORM')){

    inputs.leaf = c('N',"fqe", "Cab",  "Car", "Cs", "EWT", "LMA", "Cx","alpha")
    inputs.canopy = c("LAI", "TypeLidf", "LIDFa", "LIDFb", "hspot", "tts", "tto", "psi", "psoil",
                      "LAIu", "cd", "sd", "h", "skyl")


  } else if ((leaf.model == 'FLUSPECT-B-Cx') & (canopy.model == 'fourSAILH')){

    inputs.leaf = c('N',"fqe", "Cab",  "Car", 'Anth',"Cs",  "Cx", "EWT", "LMA","Prot","CBC","alpha")
    inputs.canopy = c("LAI", "TypeLidf", "LIDFa", "LIDFb", "hspot", "tts", "tto", "psi", "psoil")


  } else if ((leaf.model == 'FLUSPECT-B-Cx') & (canopy.model == 'INFORM')){

    inputs.leaf = c('N',"fqe", "Cab",  "Car", 'Anth',"Cs",  "Cx", "EWT", "LMA","Prot","CBC","alpha")
    inputs.canopy = c("LAI", "TypeLidf", "LIDFa", "LIDFb", "hspot", "tts", "tto", "psi", "psoil",
                      "LAIu", "cd", "sd", "h", "skyl")


  }
  # Create an empty  LUT
  LUT_db <- data.frame(matrix(ncol=length(c(inputs.leaf, inputs.canopy)), nrow=nLUT))
  colnames(LUT_db) <- c(inputs.leaf, inputs.canopy)

  # Fill LUT for each input
  for (i.input in c(inputs.leaf, inputs.canopy)) {
    min_ <- subset(LUT_ranges, inputs == i.input)$Min[1]
    max_ <- subset(LUT_ranges, inputs == i.input)$Max[1]
    set_seed <-setseed
    if (distribution == 'uniform') {

      LUT_db[[i.input]] <- stats::runif(nLUT, min=min_, max=max_)

    } else if (distribution == 'gauss') {
      mean_ <- (min_ + max_) / 2
      sd_ <- (max_ - min_) / 6  # Assuming 99.7% coverage for min/max in normal dist.
      LUT_db[[i.input]] <- stats::rnorm(nLUT, mean=mean_, sd=sd_)
    }
  }

  return(LUT_db)
}

