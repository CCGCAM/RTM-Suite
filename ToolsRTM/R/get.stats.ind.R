#' Get Statistical Scores between measurd and predicted single variable
#'
#' This function calculates statistical scores such as R-squared, RMSE, MNMB,MB, FGE and MAE.
#'
#' @param df a dataframe; the dataset.
#' @param depVar.pred A character; the predicted variable name
#' @param depVar A character; the measured variable name
#'
#' @return A data frame containing the statistical scores (R-squared, RMSE, MAE)
#' for both training and testing datasets.
#'
get.stats.ind <-function(df,depVar=NULL, depVar.pred=NULL){

  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ### Equations statistical scores
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # Function that returns Mean Absolute Error
  MAE <- function(m, o) {
    error <- m - o
    mean(abs(error), na.rm = TRUE)  # Exclude NA values
  }

  # Function that returns Root Mean Squared Error
  RMSE <- function(m, o) {
    sqrt(mean((m - o)^2, na.rm = TRUE))  # Exclude NA values
  }

  # Function that computes Mean Normalized Mean Bias (MNMB)
  MNMB <- function(observed, predicted) {
    mean_normalized_bias <- mean((predicted - observed) / observed, na.rm = TRUE)
    return(mean_normalized_bias)
  }

  # Function that computes Mean Bias (MB)
  MB <- function(observed, predicted) {
    mean_bias <- mean(predicted - observed, na.rm = TRUE)
    return(mean_bias)
  }
  # Function that computes Fractional Gross Error (FGE)
  FGE <- function(observed, predicted) {
    fractional_gross_error <- mean(abs((predicted - observed) / observed), na.rm = TRUE)
    return(fractional_gross_error)
  }




#### Skill scores for training data
r2_<-round(cor(df[,depVar.pred], df[,depVar],use='pairwise.complete.obs')^2,2)
rmse_<-round(RMSE(df[,depVar.pred], df[,depVar]),2)
mae_<-round(MAE(df[,depVar.pred], df[,depVar]),2)
mnmb_<-round(MNMB(df[,depVar.pred], df[,depVar]),2)
mb_<-round(MB(df[,depVar.pred], df[,depVar]),2)
fge_<-round(FGE(df[,depVar.pred], df[,depVar]),2)


stats<-data.frame(R2=r2_,
                  RMSE=rmse_,
                  MAE=mae_,
                  MNMB=mnmb_,
                  MB=mb_,
                  FGE=fge_)

return(stats)


}

