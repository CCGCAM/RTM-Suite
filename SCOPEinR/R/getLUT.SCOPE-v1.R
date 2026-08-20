#' Generate LUT for SCOPE (v1)
#'
#' This function generates samples of a Look-Up Table (LUT) to analyze variations
#' in plant traits required for a specific scope.
#'
#' @param inputLUT A table (LUT) containing data.
#' @param nLUT Number of samples to generate from the input LUT. Defaults to 100.
#' @param setseed Seed value to ensure reproducibility of results. Defaults to 123.
#'
#' @return Returns a data frame containing the generated LUT samples.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' inputLUT <- read.table(system.file("input", "inputs_SCOPE.csv", package = "SCOPEinR"),
#'                        header = TRUE, sep = ",")
#' getLUT.SCOPE.v1(inputLUT = inputLUT, nLUT = 100, setseed = 123)
#' }
getLUT.SCOPE.v1<-function(inputLUT, nLUT=100, setseed = 123){

  if ( is.null(setseed) | setseed ==FALSE ){
    set.seed(Sys.time())
  } else {
    set.seed(setseed)
  }

  if (is.null(nLUT)){
    message('number of varitions for each input is fixed to 100')
    nLUT = 100
  }


  if (is.null(nLUT)){
    message('number of varitions for each input is fixed to 100')
    nLUT = 100
  }

  var.list = list()
  inputs_names<- inputLUT[,'variable']


  for (i in  c(1:length(inputs_names))) {
    table.sb<-subset(inputLUT, variable == inputs_names[i])
    trait <- inputs_names[i]


    if (table.sb[,'Distribution'] == 'Uniform'){

      var.list[[trait]] <- as.numeric(stats::runif(nLUT,min = table.sb[,'lower'],max=table.sb[,'upper']))
    } else if ( table.sb[,'Distribution'] == 'Fixed'){
      var.list[[trait]] <- as.numeric(rep(table.sb[,'default'],nLUT))

    } else {
      n_casesNorm= 3 * nLUT
      var.list[[trait]] <-as.numeric(ToolsRTM::gauss_byMin_Max(n=nLUT, m=as.numeric(table.sb[,'Mean_D']),
                                                    s=as.numeric(table.sb[,'Std_D']), lwr=table.sb[,'lower'],
                                                    upr=table.sb[,'upper'], nnorm=n_casesNorm))
    }

    if(trait == 'Type') {

      var.list[[trait]] = rep(paste('C',table.sb[,'default'],sep=''),nLUT)

    }

    if(trait == 'TypeLidf') {

      # Generate a random sequence of 0s and 1s
      ramdon.value <- sample(1:2, nLUT, replace = T)

      var.list[[trait]] <- ramdon.value


    }

    ## Here important the order first TypeLidf and then LIDFa and LIDFb
    if(trait == 'LIDFa'){

      var.list[['LIDFa']]  <- ifelse(var.list[['TypeLidf']]  == 1, runif(nLUT,-1 ,1),
                                     ifelse(var.list[['TypeLidf']]  == 2, runif(nLUT,0 ,90), runif(nLUT,0 ,90 )))

    }

    if(trait == 'LIDFb'){

      var.list[['LIDFb']]  <- ifelse(var.list[['TypeLidf']]  == 1, runif(nLUT,-1 ,1),
                                     ifelse(var.list[['TypeLidf']]  == 2, 0,  runif(nLUT,-1 ,1)))

    }


    if(trait == 'startDate') {

      var.list[[trait]] = rep(paste('2018-08-01',sep=''),nLUT)

    }

    if(trait == 'endDate') {

      var.list[[trait]] = rep(paste('2018-09-01',sep=''),nLUT)

    }

  }

  ### Get correlation for Car based on Cab
  LUT <- do.call(data.frame,var.list)


  return(LUT)

}





