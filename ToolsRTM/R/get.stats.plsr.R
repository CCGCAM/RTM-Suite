#' Get Statistical Scores for PLSR Model
#'
#' This function calculates statistical scores such as R-squared, RMSE, MNMB,MB, FGE and MAE.
#' for both training and testing datasets based on predictions from a Partial Least Squares
#' Regression (PLSR) model.
#'
#' @param model The trained ML model.
#' @param k The number of components used for prediction.
#' @param train The training dataset.
#' @param test The testing dataset.
#' @param var The variable of interest.
#'
#' @return A data frame containing the statistical scores (R-squared, RMSE, MAE)
#' for both training and testing datasets.
#'
get.stats.plsr <- function(model,k,train,test,var){


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

  #################

pred.train<-predict(object = model,ncomp = k,newdata=train)
pred.test<-predict(object = model,ncomp = k,newdata=test)


#### Skill scores for training data
r2.train<-round(cor(pred.train,train[,var],use='pairwise.complete.obs')^2,2)
rmse.train<-round(RMSE(pred.train,train[,var]),2)
mae.train<-round(MAE(pred.train,train[,var]),2)
mnmb.train<-round(MNMB(pred.train,train[,var]),2)
mb.train<-round(MB(pred.train,train[,var]),2)
fge.train<-round(FGE(pred.train,train[,var]),2)

#### Skill scores for testing data
r2.test<-round(cor(pred.test,test[,var],use='pairwise.complete.obs')^2,2)
rmse.test<-round(RMSE(pred.test,test[,var]),2)
mae.test<-round(MAE(pred.test,test[,var]),2)
mnmb.test<-round(MNMB(pred.test,train[,var]),2)
mb.test<-round(MB(pred.test,train[,var]),2)
fge.test<-round(FGE(pred.test,train[,var]),2)

stats<-data.frame(R2=c(r2.train,r2.test),
                  RMSE=c(rmse.train,rmse.test),
                  MAE=c(mae.train,mae.test),
                  MNMB=c(mnmb.train,mnmb.test),
                  MB=c(mb.train,mb.test),
                  FGE=c(fge.train,fge.test))



return(stats)


}

