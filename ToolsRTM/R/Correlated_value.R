#' Generate a variable correlated with an existing one
#'
#' Draws a new variable \code{y = r*x + noise}, with noise scaled so that
#' \code{cor(x, y)} is approximately \code{r}. Negative results are clipped
#' to 0 (this assumes \code{x}/\code{y} are non-negative quantities, e.g.
#' pigment concentrations -- see \code{\link{get_distributionLUT}}'s
#' \code{DepCab} option, which uses this to correlate Car with Cab).
#'
#' @param x numeric vector. The variable to correlate against.
#' @param r numeric (-1 to 1). Target correlation coefficient between \code{x} and the result.
#'
#' @return A numeric vector the same length as \code{x}, correlated with it at approximately \code{r}.
#' @export
#'
#' @examples
#' Cab <- runif(200, 10, 80)
#' Car <- correlatedValue(x = Cab / 4, r = 0.8)
#' cor(Cab, Car)  # close to 0.8
correlatedValue = function(x, r){
  r2 = r**2
  ve = 1-r2
  SD = sqrt(ve)
  e  = rnorm(length(x), mean=0, sd=SD)
  y  = r*x + e
  ##check if some values is negative
  for (i in c(1:length(y))){
    if (y[i] <=0) {
      y[i]=0
    }
  }
  
  return(y)
}

######### Another option could be based on the relationship between Cab-Car
# constants from ANGERS03 Leaf Optical Data
#slope = 0.2234
#intercept = 0.9861
#spread = 4.6839
#grid_var=data.LUT$Cab
# car_lin = slope * grid_var+ intercept
# lower_car = slope / spread * 3 * grid_var
# upper_car = slope * spread / 3 * grid_var + 2 * intercept
# 
# car_noise=rnorm(length(grid_var), mean=3, sd=0.5)+ car_lin
