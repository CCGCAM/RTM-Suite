
#' Get LUT for radiative transfer models
#'
#' @param inputs  a table with the specific ranges, type of distribution ...
#' @param nLUT  the number of rows for the LUT
#'
#' @return a dataframe with all parameters
#' @export
#'
#'
getLUT<-function(inputs=NULL, nLUT=100, setseed = 123){

  if ( is.null(setseed)){
    set.seed(setseed)
  }
  set.seed(setseed)

  if (length(colnames(inputs)) == 10) {
    expected_columns <- c('variable', 'lower', 'upper', 'units', 'Distribution', 'Mean_D', 'Std_D', 'Dependencies', 'use.default', 'default')
    missing_columns <- setdiff(expected_columns, colnames(inputs))
  } else {
    expected_columns <- c('variable', 'lower', 'upper', 'units','model', 'Distribution', 'Mean_D', 'Std_D', 'Dependencies', 'use.default', 'default')
    missing_columns <- setdiff(expected_columns, colnames(inputs))
  }


  if (length(missing_columns) > 1) {
    message('Please provide a LUT with these names:')
    message(paste(' - ', missing_columns, collapse='\n'))
    message('Check default tables in ToolsRTM::inputs; where inputs are: inputsFlUSPECT; inputsSCOPE;, inputsSPART; inputsINFORM; inputsLiberty; inputsPROSAIL; inputsRTMs ')

    stop()

  } else if (length(missing_columns) <= 1) {
    #message('All expected columns are present in the LUT.')
  }

  if (is.null(nLUT)){
    message('number of varitions for each input is fixed to 100')
    nLUT = 100
  }

  var.list = list()
  inputs_names<- inputs[,'variable']
  for (i in  c(1:length(inputs_names))) {
    table.sb<-subset(inputs, variable == inputs_names[i])
    trait <- inputs_names[i]
    #print(inputs_names[i])
    ## for Carotenoids
    if (table.sb[,'Dependencies'] == 'Yes' ){

      var.list[[trait]] <- 'NA'
      trait_dep <-trait
    }

    if (table.sb[,'use.default'] == 1){



      if (table.sb[,'Distribution'] == 'Uniform'){

        var.list[[trait]] <- stats::runif(nLUT,min = table.sb[,'lower'],max=table.sb[,'upper'])

      } else if (table.sb[,'Distribution'] == 'Gaussian'){
        n_casesNorm= 3 * nLUT
        var.list[[trait]] <-ToolsRTM::gauss_byMin_Max(n=nLUT, m=as.numeric(table.sb[,'Mean_D']),
                                                      s=as.numeric(table.sb[,'Std_D']), lwr=table.sb[,'lower'],
                                                      upr=table.sb[,'upper'], nnorm=n_casesNorm)
      }
    } else{
      var.list[[trait]] <-table.sb[,'default']
    }
  }
  ### Get correlation for Car based on Cab
  var.list[[trait_dep]] <- ToolsRTM::correlatedValue(x=var.list[['Cab']]/4, r=.8)
  LUT <- do.call(cbind,var.list )


  return(LUT)
}
