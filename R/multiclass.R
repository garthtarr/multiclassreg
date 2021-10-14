#' Wrapper function for easy multiclass regression
#' 
#' A wrapperfunction that works with data frames, created to not ruin 
#' existing code that works directly with the multiclass_reg() function.
#'
#' @param df data frame
#' @param class_var the name of the class variable in the data frame
#' @param y_var the name of the dependent variable in the data frame
#' @param scale whether or not to scale the predictors (logical, default = TRUE)
#' @param center_y whether or not to center the response (logical, default = TRUE)
#' @param ... additional arguments to be passed on to multiclass_reg()
#' @return A list with the following components
#' \item{fit}{fitted \code{genlasso} object}
#' \item{coef}{matrix of estimated coefficients}
#' \item{lambda}{vector of regularization parameters along knots of regularization path}
#' \item{K}{number of classes}
#' \item{p}{number of predictors per class}
#' \item{n}{vector of sample sizes for each class}
#' \item{Y}{composite Y vector (stacked over the classes)}
#' \item{X}{block diagonal X matrix (each block coresponds to a class)}
#' \item{var_indicator}{matrix of variable indicators}
#' 
#' @importFrom stats coef sd
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
#' Xk <- mvtnorm::rmvnorm(n*2, mean=1:p, sigma=diag(1:p)) %>% round(2)
#' yk <- Xk %*% c(0,0,4,0,0,8,0) + rnorm(n*2) + rep(c(0,10), each = n) %>% round(2)
#' df = data.frame(y = yk, Xk, class = rep(c("A", "B"), each = n))
#' # force a missing predictor
#' df$X1[df$class == "A"] = NA
#' df$X7[df$class == "B"] = NA

#' mc = mc_reg(df)
#' mc$predictions %>% ggplot2::ggplot(ggplot2::aes(x = y, y = prediction)) +
#'   geom_abline() + 
#'   ggplot2::geom_point() +
#'   ggplot2::facet_grid(class ~ knot, labeller = label_both) +
#'   ggplot2::geom_smooth(method = "lm") +
#'   ggplot2::geom_point() + 
#'   labs(title = "Original scale before adding an intercept adjustment")

#' mc$scaled_predictions %>% ggplot2::ggplot(ggplot2::aes(x = y, y = scaled_prediction)) +
#'   geom_abline() + 
#'   ggplot2::geom_point() +
#'   ggplot2::facet_grid(class ~ knot, labeller = label_both) +
#'   ggplot2::geom_smooth(method = "lm") +
#'   ggplot2::geom_point() + 
#'   labs(title = "Scaled version")
#'   
#' mc$predictions %>% ggplot2::ggplot(ggplot2::aes(x = y, y = yhat)) +
#'  geom_abline() + 
#'  ggplot2::geom_point() +
#'  ggplot2::facet_grid(class ~ knot, labeller = label_both) +
#'  ggplot2::geom_smooth(method = "lm") +
#'  ggplot2::geom_point() + 
#'  labs(title = "Original scale with alpha added in")
#'  
#' @export
mc_reg <- function(df, 
                   class_var = "class", 
                   y_var = "y", 
                   scale = TRUE, 
                   center_y = TRUE, 
                   lambdas = NULL,
                   ...){
  
  # take a copy of the input data frame
  raw_dat = df
  
  if(center_y){
    # calculate the class means
    class_means_vector = tapply(raw_dat[[y_var]], raw_dat[[class_var]], mean, na.rm = TRUE)[raw_dat[[class_var]]]
    # center the y_var so that it has mean zero for each class
    df[y_var] = raw_dat[y_var] - class_means_vector
  } else {
    # otherwise set class_means_vector to be empty
    class_means_vector = NA
  }
  
  if(scale) {
    # don't scale the class var nor the dependent variable
    scale_cols = colnames(df)[!colnames(df)%in% c(class_var, y_var)]
    col_means = df %>% dplyr::ungroup() %>% dplyr::summarise_at(.vars = dplyr::vars(tidyselect::all_of(scale_cols)), mean, na.rm = TRUE)
    col_sds = df %>% dplyr::ungroup() %>% dplyr::summarise_at(.vars = dplyr::vars(tidyselect::all_of(scale_cols)), sd, na.rm = TRUE)
    df = df %>% dplyr::ungroup() %>% dplyr::mutate_at(.vars = dplyr::vars(tidyselect::all_of(scale_cols)), function(x) (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE))
  } else { # return the same sorts of objects but with means of zero and scale factor of 1
    scale_cols = colnames(df)[!colnames(df) %in% c(class_var, y_var)]
    col_means = df %>% dplyr::ungroup() %>% dplyr::summarise_at(.vars = dplyr::vars(tidyselect::all_of(scale_cols)), function(x) 0)
    col_sds = df %>% dplyr::ungroup() %>% dplyr::summarise_at(.vars = dplyr::vars(tidyselect::all_of(scale_cols)), function(x) 1)
    df = df %>% dplyr::ungroup()
  }
  
  lst = data_to_lists(df, class_var, y_var) 
  y = lst[[1]]
  X = lst[[2]]
  # perform the genlasso (Ines' existing code)
  out <- multiclass_reg(y, X, ...)
  # the out object has:
  # fit: a genlasso object
  # coef: a matrix of coefficients including a column at the right for 0, classes are stacked (not identified explicitly)
  # lambda: the penalty parameter values (correspond to column names in coef)
  # K: number of classes
  # p: number of variables
  # n: number of observations in each class
  # Y_list: identical to the input y
  # X_list: identical to the input X
  # Y: matrix version of Y_list
  # X: matrix version of X_list
  # var_indicator: identifying class and variable?
  
  # what was the issue with degrees of freedom?
  # these two methods seem to give the same approach
  # starts with all the same, for large lambda
  # decreases until we get to the 
  # number of variables * number of classes - 1 
  # then need to add on number of variables * number of classes
  # for the separate least squares solution
  # out$fit$df
  # out$fit$beta %>% round(5) %>% apply(2, dplyr::n_distinct)
  # I can't recall what was so bad about this?!?!
  
  ## FITTED VALUES ##
  
  # prediction for each penalty parameter (and each class)
  
  # in here, need to cbind the bls coefficients
  # like what Ines had for out$coef
  if(missing(lambdas)) {
    # ensuring distinct column names
    lambdas = unique(c(round(out$fit$lambda, 6), 0))
  } 
  out_coef = coef(out$fit, lambda = lambdas)$beta
  colnames(out_coef) = lambdas
  
  # length(unique(lambdas))
  # length(unique(colnames(out_coef)))
  
  # not sure that this will work if there are 
  # missing variables in particular classes
  unadjusted_yhat = out$X %*% out_coef #out$coef
  # when there are missing variables there's a missing column
  # in out$X and a corresponding missing row in out$coef
  # so I think this does actually work in the presence of
  # missingness....!!! does this make what I'm doing below
  # redundant?
  
  # add class_means_vector back in (intercept adjustment)
  nc = ncol(unadjusted_yhat)
  nr = nrow(unadjusted_yhat)
  if(center_y) {
    yhat = unadjusted_yhat + matrix(class_means_vector, ncol = nc, nrow = nr, byrow = FALSE)
    # store these predictions as fitted_values
    out[["fitted_values"]] = yhat
    out[["residuals"]] = matrix(raw_dat[[y_var]], ncol = nc, nrow = nr, byrow = FALSE) - yhat
  } else {
    out[["fitted_values"]] = unadjusted_yhat
    out[["residuals"]] = matrix(raw_dat[[y_var]], ncol = nc, nrow = nr, byrow = FALSE) - unadjusted_yhat
  }
  # The above works when we have missing variables in one or more classes
  
  # What we have below is a "tidy" approach which allows for (easier?)
  # identification of classes and variables, scaled and unscaled coefficients
  # and ultimately predictions too
  estimated_coef = df %>%
    dplyr::group_by(!!rlang::sym(class_var)) %>% 
    # identifies which predictors in which classes are missing
    dplyr::summarise_all(.funs = function(x) !anyNA(x)) %>% 
    tidyr::pivot_longer(cols = -tidyselect::one_of(class_var, y_var),
                        names_to = c("x_var"), 
                        values_to = "present") %>% 
    dplyr::select(- !!rlang::sym(y_var)) %>% 
    # merge in the initial scale estimate
    dplyr::left_join(tidyr::gather(col_sds, key = "x_var", value = "scale_factor"), by = "x_var") %>% 
    dplyr::filter(present) %>% 
    # bind in the coefficient matrix (assumes rows are in the correct order)
    # cbind(out$coef) %>% 
    cbind(out_coef) %>% 
    # stack so each row is a predictor-class-knot combination
    tidyr::gather(key = "lambda", 
                  value = "scaled_coef",
                  -c(class, x_var, present, scale_factor)) %>% 
    # add a knot identifier in
    # old code used group_indices
    # dplyr::mutate(knot = dplyr::group_indices(., rev(as.numeric(lambda)))) %>% 
    # new code added 4 Mar 2021 explicitly groups and ungroups
    # and use the new function cur_group_id()
    dplyr::group_by(rev(as.numeric(lambda))) %>% 
    dplyr::mutate(knot = dplyr::cur_group_id()) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(coef = scaled_coef/scale_factor) %>% 
    dplyr::group_by(knot) %>% 
    dplyr::mutate(df = dplyr::n_distinct(round(coef, 5)))
  
  # naming convention: if it doesn't have scaled in front then it means
  # unscaled original data or rescaled (i.e. transformed coefficients).
  # if it has scaled in front it means that it was built using the 
  # df data frame that scaled and/or shifted the original data
  
  pred_data = raw_dat %>% # raw dat was the original data
    # sequential observation identifier (required for matching later)
    dplyr::mutate(obs_no = 1:dplyr::n())# %>% 
    # make the class variable namified (so don't run into issues later)
    # dplyr::mutate(class = make.names(class))
  # but doing this converts it to a character, which creates issues too
  
  scaled_pred_data = df %>% # df has the scaled predictors and centered y
    # sequential observation identifier (required for matching later)
    dplyr::mutate(obs_no = 1:dplyr::n())# %>% 
    # make the class variable namified (so don't run into issues later)
    # dplyr::mutate(class = make.names(class))
  # but doing this converts it to a character, which creates issues too
  
  # extract a data frame with the y data
  pred_y_data = pred_data[,c("obs_no", y_var, class_var)]
  scaled_pred_y_data = scaled_pred_data[,c("obs_no", y_var, class_var)]
  
  # define the x testing data frame (removing the y variable)
  pred_x_data = pred_data %>% 
    dplyr::select(- !!rlang::sym(y_var)) %>% 
    # convert the x data frame into long form test_x_data_df
    tidyr::gather(key = x_var,
                  value = x_value, -class, -obs_no)
  
  scaled_pred_x_data = scaled_pred_data %>% 
    dplyr::select(- !!rlang::sym(y_var)) %>% 
    # convert the x data frame into long form test_x_data_df
    tidyr::gather(key = x_var,
                  value = x_value, -class, -obs_no)
  
  # calculate the residuals by joining the x data in
  # with the coefficient date frame, this repeats the
  # estiamted coefficients down the data frame for the 
  # number of observations in each variable and class
  preds = dplyr::left_join(estimated_coef, 
                           pred_x_data, 
                           by = c("class", "x_var")) %>% 
    dplyr::group_by(obs_no, class, knot, df) %>% 
    dplyr::summarise(prediction = sum(coef * x_value)) %>% 
    dplyr::left_join(pred_y_data, by = c("obs_no", "class"))
  
  scaled_preds = dplyr::left_join(estimated_coef, 
                                  scaled_pred_x_data, 
                                  by = c("class", "x_var")) %>% 
    dplyr::group_by(obs_no, class, knot) %>% 
    dplyr::summarise(scaled_prediction = sum(scaled_coef * x_value)) %>% 
    dplyr::left_join(scaled_pred_y_data, 
                     by = c("obs_no", "class"))
  
  alpha = preds %>% dplyr::group_by(class, knot) %>% 
    dplyr::summarise(class_means = mean(!!rlang::sym(y_var)),
                     scaled_pred_means = mean(prediction)) %>% 
    dplyr::mutate(alpha = class_means - scaled_pred_means) %>% 
    dplyr::select(-class_means, -scaled_pred_means)
  
  preds = preds %>% dplyr::left_join(alpha, by = c("class", "knot")) %>% 
    dplyr::mutate(yhat = prediction + alpha) 
  
  # probably should just remove the "prediction" from preds object
  # but will leave it in temporarily
  
  gamma = 0.5
  Kp = max(preds$df)
  bic = preds %>% dplyr::group_by(knot, df) %>% 
    dplyr::summarise(
      sigma2 = mean((y-yhat)^2),
      df = mean(df),
      n = dplyr::n(),
      aic = log(sigma2) + df * 2 / n,
      bic1 = log(sigma2) + df * log(n) / n,
      bic2 = sigma2 + df * log(n) / n,
      # update the ebic
      # something strange going on with the two bic
      ebic = log(sigma2) + df * (log(n) + 2*gamma*log(Kp))/n
    )
  
  # now need to make the alpha coefficients accessible to 
  # external prediction
  
  out[["coef_df"]] = estimated_coef
  out[["raw_data"]] = raw_dat
  out[["class_var"]] = class_var
  out[["y_var"]] = y_var
  out[["scale"]] = scale
  out[["center_y"]] = center_y
  out[["class_means_vector"]] = class_means_vector
  out[["alpha"]] = alpha
  out[["predictions"]] = preds
  out[["scaled_predictions"]] = scaled_preds
  out[["lambdas"]] = lambdas
  out[["bic"]] = bic
  
  # define a class for plotting
  class(out) = c("mcreg")
  
  return(out)
}

#' Predict method for mcreg objects
#' 
#' Predict given an input.
#'
#' @param obj an multiclass regression object (object of class mcreg)
#' @param class_var the name of the class variable in the data frame
#' @param y_var the name of the dependent variable in the data frame
#' @param ... additional arguments to be passed on to multiclass_reg()
#' @return A list with the following components
#' @examples
#' p = 7
#' k = 2
#' n = 20
#' beta = c(1,2)
#' set.seed(1)
#' X = list(matrix(rnorm(p*n), ncol = p), matrix(rnorm(p*n), ncol = p))
#' Y = list(rnorm(n), rnorm(n))
#' df = lists_to_data(Y, X)
#' multiclass_reg(df)
#' @export
predict.mcreg <- function(obj, 
                          class_var = "class", 
                          x_vars = NULL, 
                          coef_type = c("rescaled","scaled"),
                          plot_type = c("knot","lambda"),
                          log_x = NULL){
  
}

#' Plot coefficient paths for multiclass regression objects
#' 
#' Plot the coefficient paths.
#'
#' @param obj an multiclass regression object (object of class mcreg)
#' @param class_var the name of the class variable in the data frame
#' @param y_var the name of the dependent variable in the data frame
#' @param ... additional arguments to be passed on to multiclass_reg()
#' @return A list with the following components
#' \item{fit}{fitted \code{genlasso} object}
#' \item{coef}{matrix of estimated coefficients}
#' \item{lambda}{vector of regularization parameters along knots of regularization path}
#' \item{K}{number of classes}
#' \item{p}{number of predictors per class}
#' \item{n}{vector of sample sizes for each class}
#' \item{Y}{composite Y vector (stacked over the classes)}
#' \item{X}{block diagonal X matrix (each block coresponds to a class)}
#' \item{var_indicator}{matrix of variable indicators}
#' @examples
#' p = 7
#' k = 2
#' n = 20
#' beta = c(1,2)
#' set.seed(1)
#' X = list(matrix(rnorm(p*n), ncol = p), matrix(rnorm(p*n), ncol = p))
#' Y = list(rnorm(n), rnorm(n))
#' df = lists_to_data(Y, X)
#' mc_reg(df)
#' @export
plot.mcreg <- function(obj, 
                       class_var = "class", 
                       x_vars = NULL, 
                       coef_type = c("rescaled","scaled"),
                       plot_type = c("knot","lambda"),
                       log_x = NULL){
  
  y_var = switch(coef_type[1], 
                 "rescaled" = "coef",
                 "scaled" = "scaled_coef")
  
  y_lab = switch(coef_type[1], 
                 "rescaled" = "Rescaled coefficients",
                 "scaled" = "Standardised coefficients")
  
  x_var = switch(plot_type[1],
                 "knot" = "knot",
                 "lambda" = "lambda")
  
  x_lab = switch(plot_type[1],
                 "knot" = "Knot",
                 "lambda" = "Lambda")
  
  obj$coef_df$knot = as.numeric(obj$coef_df$knot)
  obj$coef_df$lambda = as.numeric(obj$coef_df$lambda)
  
  if(missing(log_x)) {
    if(plot_type[1] == "lambda") {
      log_x = TRUE
    } else {
      log_x = FALSE
    }
  }
  
  if(missing(x_vars)) {
    x_vars = unique(obj$coef_df$x_var)
  } else if (!all(x_vars %in% unique(obj$coef_df$x_var))) {
    stop("Invalid 'x_vars' value")
  }
  
  p = obj$coef_df %>% 
    dplyr::filter(x_var %in% x_vars) %>% 
    ggplot2::ggplot(ggplot2::aes_string(x = x_var, y = y_var, colour = class_var, group = class_var)) + 
    ggplot2::geom_hline(yintercept=0, alpha = 0.5) +
    # geom_vline(xintercept = cvfit$l1opt.knot, alpha = 0.5)+
    ggplot2::geom_line() + 
    ggplot2::facet_wrap(~x_var, scales = "free_y") + 
    ggplot2::scale_color_brewer(palette = "Set1") + 
    ggplot2::labs(x = x_lab, y = y_lab, colour = "Class") + 
    ggplot2::theme_classic()
  
  if(x_var == "lambda"){
    p = p + ggplot2::scale_x_continuous(trans=scales::pseudo_log_trans(base = 10))
  }
  
  return(p)
  
}

#' Fit the multiclass regression estimator
#'
#' @param Y list of size K. Its kth entry is a vector of size n_k
#' @param X list of size K. Its kth entry is a matrix of size n_kxp
#' @param weighting logical. If TRUE: use weighted regularization parameter
#' @param weights provide vector of weights
#' @param wtype either "size" (default, weighting based on group samle sizes), "init" (weighting based on initial estimator), "combi" (weighting based on group sample sizes and initial estimator)
#' @param initestim pK-dimensional vector with estimated coefficients of initial estimator
#' @param normalize logical. Whether or not to normalize the weights (currently divides by maximum).
#' @return A list with the following components
#' \item{fit}{fitted \code{genlasso} object}
#' \item{coef}{matrix of estimated coefficients}
#' \item{lambda}{vector of regularization parameters along knots of regularization path}
#' \item{K}{number of classes}
#' \item{p}{number of predictors per class}
#' \item{n}{vector of sample sizes for each class}
#' \item{Y}{composite Y vector (stacked over the classes)}
#' \item{X}{block diagonal X matrix (each block coresponds to a class)}
#' \item{var_indicator}{matrix of variable indicators}
#' @examples
#' p = 7
#' k = 2
#' n = 20
#' beta = c(1,2)
#' set.seed(1)
#' X = list(matrix(rnorm(p*n), ncol = p), matrix(rnorm(p*n), ncol = p))
#' Y = list(rnorm(n), rnorm(n))
#' # no missing data
#' res = multiclass_reg(Y, X)
#' @export
multiclass_reg <- function(
  Y, X,
  weighting = FALSE,
  weights = NULL,
  wtype = "size", 
  initestim = NULL, 
  normalize = TRUE,
  x_names = NULL,
  lambdas = NULL # not implemented
  ){
  
  # Depends on Matrix, genlasso package, Uses inputs functions
  
  # Initial checks
  if(missing(x_names)){
    x_names = paste0("V",1:ncol(X[[1]]))
  }
  
  if(!class(Y)=="list"){
    stop("Please provide a list object for Y")
  }
  
  if(!class(X)=="list"){
    stop("Please provide a list object for X")
  }
  
  if(!length(Y)==length(X)){
    stop("The number of classes in Y and X do not match")
  }
  
  # if(length(unique(unlist(lapply(X, function(U){ncol(U)}))))!=1){
  #   stop("The number of predictors do not match across classes")
  # }
  
  if(wtype == "init" & is.null(initestim)){
    stop("Provide an initial estimator under initestim if wtype is set to init")
  }
  
  if(wtype == "combi" & is.null(initestim)){
    stop("Provide an initial estimator under initestim if wtype is set to combi")
  }
  
  # Start code
  # Dimensions
  K <- length(Y) # number of classes
  pseq <- sapply(X, ncol) # number of predictors in each class
  p <- max(pseq)
  n <-  sapply(Y, length)
  
  # Get inputs
  get.input <- inputs(Y = Y, X = X, 
                      weighting = weighting,
                      weights = weights,
                      wtype = wtype, 
                      buildD = TRUE,
                      normalize = normalize)
  
  # Multi-class fit
  # want to optionally specify lambda sequence here
  # not sure that we can?!?
  fit <- genlasso::genlasso(y=get.input$Y,
                            X=get.input$X,
                            D=get.input$D)
  
  # Store results
  bhats <- cbind(fit$beta, fit$bls)
  colnames(bhats) <- c(colnames(fit$beta), 0)
  rownames(bhats) = colnames(get.input$X)
  
  # bhats <- unique(bhats, MARGIN=2) # Estimated coefficients
  lambdas <- as.numeric(colnames(bhats)) # regularization parameters
  
  # THINK ABOUT MOVING ALL THIS TO THE mc_reg() FUNCTION
  
  # preferable to have a different format for plotting coefficient paths
  # coef_paths = data.frame(bhats)
  # colnames(coef_paths) = paste0("Knot_",1:length(lambdas))
  # # coef_paths = tibble::rownames_to_column(data.frame(coef_paths), var = "var_name")
  # var_names = unlist(lapply(X, FUN = function(x) colnames(x)[apply(x,2,FUN = function(y) !all(is.na(y)))]))
  # class_p = unlist(lapply(X, FUN = function(x) length(colnames(x)[apply(x,2,FUN = function(y) !all(is.na(y)))])))
  # class_names = rep(names(class_p),class_p)
  # coef_paths$var_names = var_names
  # coef_paths$class_names = class_names
  # coef_long = tidyr::gather(coef_paths, key = "knot", value = "coef", -class_names, -var_names)
  # coef_long$lambda = rep(c(fit$lambda,0),each = nrow(coef_paths))
  
  out <- list("fit" = fit,
              "coef" = bhats,
              "lambda" = lambdas,
              "K" = K, "p" = p, "n" = n,
              "Y_list" = Y,
              "X_list" = X,
              "Y" = get.input$Y,
              "X" = get.input$X,
              "var_indicator" = get.input$var_indicator#,
              # "coef_paths" = coef_paths,
              # "coef_long" = coef_long
  )
}
