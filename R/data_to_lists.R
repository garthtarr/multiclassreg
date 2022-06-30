#' Data frame to lists
#'
#' @description Takes a data frame and return the lists required for \code{muliclass_reg()} to work
#' @param data data frame with (at least) all the variables in the formula as 
#'   well as the class identifier variable
#' @param class_var the name of the class variable in the data frame
#' @param y_var the name of the dependent variable in the data frame
#' @return A list with the following components
#' \item{Y}{list of size \eqn{K}. Its kth entry is a vector of size \eqn{n_k}}
#' \item{X}{list of size \eqn{K}. Its kth entry is a matrix of size \eqn{n_{k} \times p}}
#' 
#' @examples
#' p = 7
#' k = 2
#' n = 20
#' beta = c(1,2)
#' set.seed(1)
#' X = list(matrix(rnorm(p*n), ncol = p), matrix(rnorm(p*n), ncol = p))
#' Y = list(rnorm(n), rnorm(n))
#' df = lists_to_data(Y, X)
#' lst = data_to_lists(df)
#' 
#' @export
data_to_lists <- function(
  data, 
  class_var = "class", 
  y_var = "y"
  ){
  
  y = split(data[,y_var, drop = FALSE], f = data[,class_var])
  y = lapply(y, as.matrix)
  
  X = split(data[,!names(data) %in% c(y_var, class_var)], f = data[,class_var])
  X = lapply(X, as.matrix)

  return(list(y = y, X = X))
}

#' List to data frame
#'
#' @description Take lists and convert them to data frame
#' @param Y list of response vectors
#' @param X list of predictor matrices
#' @return A data frame with all the X variables, 
#' as well as the Y variable and a new column identifying 
#' the class.
#' 
#' @examples
#' p = 7
#' k = 2
#' n = 20
#' beta = c(1,2)
#' set.seed(1)
#' X = list(matrix(rnorm(p*n), ncol = p), matrix(rnorm(p*n), ncol = p))
#' Y = list(rnorm(n), rnorm(n))
#' df = lists_to_data(Y, X)
#' lst = data_to_lists(df)
#'  
#' @export
lists_to_data = function(Y, X){
  if(length(X)!=length(Y)){
    stop("Number of classes in X does not equal the number of classes in Y.")
  }
  Xdf = dplyr::bind_rows(lapply(X, as.data.frame), .id = "class")
  Ydf = dplyr::bind_rows(lapply(Y, as.data.frame))
  colnames(Ydf) = "y"
  df = dplyr::bind_cols(Xdf, Ydf)
  return(df)
}
