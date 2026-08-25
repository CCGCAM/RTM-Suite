# ------------------------------------------------------------
# Function: compute_stats()
# Returns a compact table of model performance metrics
# ------------------------------------------------------------

#' Compute statistical metrics for model performance
#'
#' This function evaluates training and testing performance
#' using common statistical accuracy metrics such as RMSE,
#' MAE, MB, MNMB, FGE, and R².
#'
#' @param pred.train Numeric vector of model predictions (training set)
#' @param train.obs  Numeric vector of observed values (training set)
#' @param pred.test  Numeric vector of model predictions (testing set)
#' @param test.obs   Numeric vector of observed values (testing set)
#'
#' @return A data.frame with Train/Test RMSE, MAE, MB, MNMB, FGE, and R².
#' @export
#'
#' @examples
#' \dontrun{
#' stats <- compute_stats(pred.train, train.obs, pred.test, test.obs)
#' }
compute_stats <- function(pred.train, train.obs,
                          pred.test,  test.obs) {

  MAE <- function(m, o) {
    error <- m - o
    mean(abs(error), na.rm = TRUE)
  }
  RMSE <- function(m, o) {
    sqrt(mean((m - o)^2, na.rm = TRUE))
  }
  MNMB <- function(observed, predicted) {
    mean_normalized_bias <- mean((predicted - observed)/observed,
                                 na.rm = TRUE)
    return(mean_normalized_bias)
  }
  MB <- function(observed, predicted) {
    mean_bias <- mean(predicted - observed, na.rm = TRUE)
    return(mean_bias)
  }
  FGE <- function(observed, predicted) {
    fractional_gross_error <- mean(abs((predicted - observed)/observed),
                                   na.rm = TRUE)
    return(fractional_gross_error)
  }


  stats <- data.frame(
    Dataset = c("Train", "Test"),

    RMSE = c(
      round(RMSE(pred.train, train.obs), 3),
      round(RMSE(pred.test,  test.obs), 3)
    ),

    MAE = c(
      round(MAE(pred.train, train.obs), 3),
      round(MAE(pred.test,  test.obs), 3)
    ),

    MNMB = c(
      round(MNMB(pred.train, train.obs), 3),
      round(MNMB(pred.test,  test.obs), 3)
    ),

    MB = c(
      round(MB(pred.train, train.obs), 3),
      round(MB(pred.test,  test.obs), 3)
    ),

    FGE = c(
      round(FGE(pred.train, train.obs), 3),
      round(FGE(pred.test,  test.obs), 3)
    ),

    R2 = c(
      round(cor(pred.train, train.obs, use = "pairwise.complete.obs")^2, 3),
      round(cor(pred.test,  test.obs,  use = "pairwise.complete.obs")^2, 3)
    )
  )

  return(stats)
}
