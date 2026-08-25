

#' get random plant traits using a uniform distribution
#'
#' @param RTmodel RT model: options are: 'PROSPECT-PRO';'PROSPECT-D'; 'INFORM'; and 'fourSAILH'
#' @param nLUT nummber of combinations
#' @param random.input  the seed number for the  uniform distribution
#'
#' @return a LUT 
#' @export
#'
#' 
#' 
get_data <- function(RTmodel='PROSPECT-D', nLUT=100, random.input =100) {
  
  if (is.null(nLUT)){
    message('number of varitions for each input is fixed to 100')
    nLUT = 100
  }
  
  if (is.null(random.input)){
    message('number of varitions for each input is fixed to 100')
    random.input = 1234
    set.seed(random.input)
  } else {
    set.seed(random.input)
  }

  
  if (RTmodel == "PROSPECT-PRO") {
    # Create lookup table for prospect model
    lut <- data.frame(
      N = runif(nLUT, min = 1, max = 3),
      Cab = runif(nLUT, min = 0, max = 100),
      Car = runif(nLUT, min = 0, max = 40),
      Anth = runif(nLUT, min = 0, max = 7.5),
      Cbrown = runif(nLUT, min = 0, max = 1),
      EWT = runif(nLUT, min = 0.001, max = 0.04),
      LMA = runif(nLUT, min = 0.001, max = 0.04),
      alpha = runif(nLUT, min = 10, max = 50),
      Prot = runif(nLUT, min = 0.001, max = 0.02),
      CBC = runif(nLUT, min = 0.001, max = 0.02)
    )
  } else if (RTmodel == "PROSPECT-D") {
      # Create lookup table for prospect model
      lut <- data.frame(
        N = runif(nLUT, min = 1, max = 3),
        Cab = runif(nLUT, min = 0, max = 100),
        Car = runif(nLUT, min = 0, max = 40),
        Anth = runif(nLUT, min = 0, max = 7.5),
        Cbrown = runif(nLUT, min = 0, max = 1),
        EWT = runif(nLUT, min = 0.001, max = 0.04),
        LMA = runif(nLUT, min = 0.001, max = 0.04),
        alpha = runif(nLUT, min = 10, max = 50)
      )

    } else if (RTmodel == "Liberty") {
      # Create lookup table for liberty model
      lut <- data.frame(
        cell.d = runif(nLUT, min = 20, max = 200),
        inter.c = runif(nLUT, min = 0.01, max = 0.1),
        baseline.abs = runif(nLUT, min = 0.0004, max = 0.0006),
        leaf.thick = runif(nLUT, min = 1, max = 10),
        albino.abs = runif(nLUT, min = 0, max = 4),
        Cab = runif(nLUT, min = 0, max = 100),
        EWT = runif(nLUT, min = 0.001, max = 0.04),
        lign.cell = runif(nLUT, min = 10, max = 80),
        Nitrogen = runif(nLUT, min = 0.3, max = 2)
 
      )

  } else if (RTmodel == "fourSAILH") {
    # Create lookup table for fourSAILH model
    lut <- data.frame(
      
      LAI  = runif(nLUT, min = 0, max = 9),
      TypeLidf = 2,
      LIDFa = runif(nLUT, min = 0, max = 90),
      LIDFb = 0,
      hspot  = runif(nLUT, min = 0, max = 1),
      tts = runif(nLUT, min = 0, max = 30),
      tto = runif(nLUT, min = 0, max = 55),
      psi = runif(nLUT, min = 0, max = 180),
      psoil = runif(nLUT, min = 0, max = 1)
     
    )
  } else if (RTmodel == "INFORM") {
    # Create lookup table for INFORM model
    lut <- data.frame(
      LAI  = runif(nLUT, min = 0, max = 9),
      TypeLidf = 2,
      LIDFa = runif(nLUT, min = 0, max = 90),
      LIDFb = 0,
      hspot  = runif(nLUT, min = 0, max = 1),
      tts = runif(nLUT, min = 0, max = 30),
      tto = runif(nLUT, min = 0, max = 55),
      psi = runif(nLUT, min = 0, max = 180),
      psoil = runif(nLUT, min = 0, max = 1),
      LAIu = runif(nLUT, min = 0.0, max = 2),
      cd = runif(nLUT, min = 0, max = 10),
      sd = runif(nLUT, min = 0, max = 1200),
      h = runif(nLUT, min = 1, max = 40),
      skyl = runif(nLUT, min = 0, max = 0.2)
    )
 
  } else {
    # Throw an error if an invalid model is specified
    stop("Invalid model specified.")
  }
 
  return(lut)
}
