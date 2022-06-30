# To do (future?)



#' Predict method for mcreg objects
#' 
#' Predict given an input.
#'
#' @param obj an multiclass regression object (object of class mcreg)
#' @param class_var the name of the class variable in the data frame
#' @param y_var the name of the dependent variable in the data frame
#' @param ... additional arguments to be passed on to multiclass_reg()
#' @return A list with the following components
#' @export
predict.mcreg <- function(obj, 
                          class_var = "class", 
                          x_vars = NULL, 
                          coef_type = c("rescaled","scaled"),
                          plot_type = c("knot","lambda"),
                          log_x = NULL){
  
}