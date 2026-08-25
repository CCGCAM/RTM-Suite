#################################
###########

#' Calculate Variance Inflation Factor (VIF)
#'
#' This function calculates the Variance Inflation Factor (VIF) for the predictor variables
#' in a linear regression model to assess multicollinearity.
#'
#' @param in_frame A data frame containing the predictor variables. The dependent variable
#'                 should not be included in this frame.
#' @param thresh A numeric value indicating the threshold for VIF; predictors with VIF
#'               greater than this threshold will be flagged as having multicollinearity.
#'               Default is 10.
#' @param trace A boolean value; if TRUE, the function will print information about
#'              the VIF calculations and any predictors that exceed the threshold.
#' @param ... Additional arguments passed to other methods (not used in this function).
#'
#' @return A data frame containing the VIF values for each predictor variable.
#' @export
#'
#' @examples
#' # Example data frame
#' df <- data.frame(x1 = rnorm(100), x2 = rnorm(100), x3 = rnorm(100))
#' df$x2 <- df$x1 + rnorm(100, sd = 0.1)  # Introduce multicollinearity
#' vif_results <- getVIF(in_frame = df, thresh = 5, trace = TRUE)
#' print(vif_results)
#' # The following VIF function were extracted from 
#' #https://beckmw.wordpress.com/2013/02/05/collinearity-and-stepwise-vif-selection/
#' 
getVIF<-function(in_frame,thresh=10,trace=T,...){

  if(!inherits(in_frame, 'data.frame')) in_frame<-data.frame(in_frame)
  
  #get initial vif value for all comparisons of variables
  vif_init<-NULL
  var_names <- names(in_frame)
  for(val in var_names){
    regressors <- var_names[-which(var_names == val)]
    form <- paste(regressors, collapse = '+')
    form_in <- formula(paste(val, '~', form))
    vif_init<-rbind(vif_init, c(val, fmsb::VIF(lm(form_in, data = in_frame, ...))))
  }
  vif_max<-max(as.numeric(vif_init[,2]))
  
  if(vif_max < thresh){
    if(trace==T){ #print output of each iteration
      prmatrix(vif_init,collab=c('var','vif'),rowlab=rep('',nrow(vif_init)),quote=F)
      cat('\n')
      cat(paste('All variables have VIF < ', thresh,', max VIF ',round(vif_max,2), sep=''),'\n\n')
    }
    return(var_names)
  }
  else{
    
    in_dat<-in_frame
    
    #backwards selection of explanatory variables, stops when all VIF values are below 'thresh'
    while(vif_max >= thresh){
      
      vif_vals<-NULL
      var_names <- names(in_dat)
      
      for(val in var_names){
        regressors <- var_names[-which(var_names == val)]
        form <- paste(regressors, collapse = '+')
        form_in <- formula(paste(val, '~', form))
        vif_add<-fmsb::VIF(lm(form_in, data = in_dat, ...))
        vif_vals<-rbind(vif_vals,c(val,vif_add))
      }
      max_row<-which(vif_vals[,2] == max(as.numeric(vif_vals[,2])))[1]
      
      vif_max<-as.numeric(vif_vals[max_row,2])
      
      if(vif_max<thresh) break
      
      if(trace==T){ #print output of each iteration
        prmatrix(vif_vals,collab=c('var','vif'),rowlab=rep('',nrow(vif_vals)),quote=F)
        cat('\n')
        cat('removed: ',vif_vals[max_row,1],vif_max,'\n\n')
        flush.console()
      }
      
      in_dat<-in_dat[,!names(in_dat) %in% vif_vals[max_row,1]]
      
    }
    
    return(names(in_dat))
    
  }
  
  
}
