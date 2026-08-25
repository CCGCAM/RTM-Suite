#' This function plots predictions from a machine learning (ML) model
#' on both training and testing datasets.
#'
#' @param model The trained ML model.
#' @param data.train The training dataset, must include the variable of interest.
#' @param data.test The testing dataset, must include the variable of interest.
#' @param var The variable of interest to plot (as a string).
#'
#' @return A ggplot object containing the plot of predictions.
#' @export
#'
#' @examples
#' \dontrun{
#' # Example usage
#' library(ggplot2)
#' # Assuming 'model', 'data.train', and 'data.test' are defined and 'var' is the target variable
#' plot <- get.plot.ML(model, data.train, data.test, var = "target_variable")
#' print(plot)
#' 
#' }
get.plot.ML<-function(model,data.train, data.test, var){

  pred.train<-predict(object = model,newdata=data.train)
  pred.test<-predict(object = model,newdata=data.test)
  stats<-get.stats(model,data.train,data.test,var)


  preds <-cbind (data.test,Preds= pred.test)

  mylabel.r = bquote(bold(r)^2 == .(format(stats[1,1], digits = 3)))
  mylabel.rmse = bquote(bold(rmse) == .(format(stats[1,2], digits = 3)))
  mylabel.r.test = bquote(bold(r)^2 == .(format(stats[2,1], digits = 3)))
  mylabel.rmse.test = bquote(bold(rmse) == .(format(stats[2,2], digits = 3)))


  axis_y<-bquote(bold(.(var)['measured'])) # axis x
  axis_x<-bquote(bold(.(var)['predicted']))# axis y

  # Create plot
  plot.ind <- ggplot(preds, aes_string(x = var, y = 'Preds')) +
    geom_point(size = 1, alpha = 0.5) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
    labs(title = paste("ML model: R^2 =", round(stats[2, 1], digits = 3), ", RMSE =", round(stats[2, 2], digits = 3)),
         x = axis_x, y = axis_y) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
          panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = 'white', color = 'white'),
          legend.key = element_rect(fill = "white", color = "white"),
          axis.title = element_text(face = "bold", size = 12),
          axis.text.y = element_text(hjust = 0.5, size = 12, face = "bold"),
          axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 12, face = "bold"),
          legend.title = element_blank()) +   theme_bw()
  print(plot.ind)

}
