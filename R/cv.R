#' New function to compute optimal value of regularization parameter based on K-fold cross-validation
#' @param fit fitted object returned from mc_reg()
#' @param Kfold K-fold cross-validation
#' @param nweight vector of length K indicating weights for MSFE measure
#' @param asym flag: asymmetric forecast error measure or not
#' @param type MSFE or MAFE
#' @param l1length length sparsity grid l1 penalty
#' @return A list with the following components
#' \item{lambda_opt}{Selected value of the regularization parameter via BIC}
#' \item{bhat_opt}{Estimated coefficients selected via BIC}
#' @examples
#' p <- 7
#' k <- 2
#' n <- 500
#' beta <- c(1, 2)
#' set.seed(1)
#' X <- list(matrix(rnorm(p * n), ncol = p), matrix(rnorm(p * n), ncol = p))
#' Y <- list(rnorm(n), rnorm(n))
#' df = lists_to_data(Y, X)
#' res = mc_reg(df)
#' cv_res <- mc_cv(res)
#' 
#' p = 7
#' k = 2
#' n = 100
#' beta = c(1,2)
#' set.seed(2)
#' X = list(matrix(rnorm(p*n), ncol = p), matrix(rnorm(p*n), ncol = p))
#' Y = list(rnorm(n), rnorm(n))
#' df = lists_to_data(Y, X)
#' Xk <- mvtnorm::rmvnorm(n*2, mean=1:p, sigma=diag(1:p)) %>% round(2)
#' yk <- Xk %*% c(0,0,4,0,0,8,0) + rnorm(n*2) + rep(c(0,100), each = n) %>% round(2)
#' df = data.frame(y = yk, Xk, class = rep(c("A", "B"), each = n))
#' res = mc_reg(df)
#' cv_res <- mc_cv(res)
#' 
#' # different cv runs have different lamda values
#' # is this a problem?
#' cv_res$cv_fits[[1]]$lambda
#' cv_res$cv_fits[[2]]$lambda
#' # should fix the lambda grid across all lambda runs
#' # with fixed lambda values, we need to force that
#' # using either the full original data 
#' # Note this is "fixed" because we don't use the 
#' # default lambda grid when extracting coefficients
#' 
#' 
#' p1 = cv_res$class_knot_summary %>%
#'   ggplot2::ggplot(ggplot2::aes(x = knot, colour = class)) +
#'   ggplot2::geom_point(ggplot2::aes(y = MSFE)) +
#'   ggplot2::geom_line(ggplot2::aes(y = MSFE))
#'
#' p2 = cv_res$class_knot_summary %>%
#'   ggplot2::ggplot(ggplot2::aes(x = knot, colour = class)) +
#'   ggplot2::geom_point(ggplot2::aes(y = MAFE)) +
#'   ggplot2::geom_line(ggplot2::aes(y = MAFE))
#'   
#' p3 = cv_res$knot_summary %>%
#'  ggplot2::ggplot(ggplot2::aes(x = knot, group = "")) +
#'  ggplot2::geom_point(ggplot2::aes(y = MSFE), colour = "blue") +
#'  ggplot2::geom_point(ggplot2::aes(y = MAFE), colour = "red") +
#'  ggplot2::geom_line(ggplot2::aes(y = MSFE), colour = "blue") +
#'  ggplot2::geom_line(ggplot2::aes(y = MAFE), colour = "red")
#'  
#'  gridExtra::grid.arrange(p1, p2, p3)
#' @export
mc_cv <- function(fit, # mcreg object
                  Kfold = 5,
                  nweight = NULL,
                  weighting = TRUE,
                  wtype = "size",
                  type = "MSFE",
                  l1length = 100,
                  normalize = TRUE) {
  
  # can run existing CV on this object
  # cv_old = multiclass_cv(fit, Kfold = Kfold,
  #                        nweight = nweight, weighting = weighting,
  #                        wtype = wtype, type = type,
  #                        l1length = l1length, normalize = normalize)
  
  # new CV approach
  
  # assign observations to fold groups
  group = fit$raw_data[[fit$class_var]]
  fold_id = caret::createFolds(y = factor(group), 
                               k = Kfold, 
                               list = FALSE)
  
  # pooled full model: drop any variables with missingness because 
  # can't predict with the missing observations
  complete_raw_dat = fit$raw_data[,!sapply(fit$raw_data, anyNA)]
  complete_names = names(complete_raw_dat)
  x_vars = paste(complete_names[!complete_names %in% c(fit$y_var)], collapse = "+")
  # x_vars = paste(complete_names[!complete_names %in% c(fit$y_var,fit$class_var)], collapse = "+")
  pooled_formula = formula(paste(fit$y_var, "~", x_vars))
  
  # initialise result storage objects
  pooled_cv_fits = cv_fits = train = test = list()
  
  # run the mc_reg() function over the Kfolds and
  # capture the results
  for(fold in seq_len(Kfold)){
    train[[fold]] = fit$raw_data[fold_id != fold,]
    test[[fold]] = fit$raw_data[fold_id == fold,]
    # note the problem here, there's different scaling for each fold
    # kind of want to hold the scaling constant across the folds
    
    # to do: force lambda grid in the cv_fits,
    # this is passed in now, but not implemented
    # can't see how to feed it into genlasso
    cv_fits[[fold]] = mc_reg(train[[fold]], 
                             class_var = fit$class_var,
                             y_var = fit$y_var,
                             scale = fit$scale,
                             center_y = fit$center_y,
                             lambdas = fit$lambdas)
    
    pooled_cv_fits[[fold]] = lm(pooled_formula, data = train[[fold]])
  }
  
  # OK now we need to evaluate performance for each run
  # we need to do this on the standardised scale
  # as well as the original scale
  
  # the coef_df table has the coefficients on the standardised
  # scale as well as the rescaled coefficients
  
  # head(fit$coef_df)
  
  # generate standardised and rescaled predictions for the testing data set
  # pull out the X matrix for each class in the test set:
  # data_to_lists(fit$coef_df, class_var = )
  # data_to_lists(test[[1]], class_var = fit$class_var, y_var = fit$y_var)$X
  
  pooled_forecast_errors = forecast_errors = list()
  
  # think about converting to multicore loop
  for(fold in seq_len(Kfold)){
    test_data = test[[fold]]
    # sequential observation identifier (required for matching later)
    test_data$obs_no = 1:nrow(test_data)
    
    # make the class variable namified (so don't run into issues later)
    # test_data$class = make.names(test_data$class)
    # doing this means we end up with character instead of factor
    
    # extract a data frame with the y data
    test_y_data = test_data[,c("obs_no", fit$y_var, fit$class_var)]
    # define the x testing data frame (removing the y variable)
    test_x_data = test_data
    test_x_data[[fit$y_var]] = NULL
    # convert the x data frame into long form
    test_x_data_df = tidyr::gather(test_x_data, key = x_var,
                                   value = x_value, -class, -obs_no)
    # calculate the forecast errors by joining the x data in
    # with the coefficient date frame, this repeats the
    # estiamted coefficients down the data frame for the 
    # number of observations in each variable and class
    forecast_errors[[fold]] = dplyr::left_join(cv_fits[[fold]]$coef_df, 
                                            test_x_data_df, 
                                            by = c("class","x_var")) %>% 
      dplyr::group_by(obs_no, class, knot) %>% 
      dplyr::summarise(prediction = sum(coef * x_value)) %>% 
      dplyr::ungroup() %>% 
      dplyr::left_join(test_y_data, by = c("obs_no", "class")) %>% 
      dplyr::left_join(cv_fits[[fold]]$alpha, by = c("class", "knot")) %>% 
      dplyr::mutate(yhat = prediction + alpha) %>% 
      dplyr::mutate(fe = !!rlang::sym(fit$y_var) - yhat)
    
    classes = test_data[[fit$class_var]]
    errors = predict(pooled_cv_fits[[fold]], newdata = test_data) - test_data[[fit$y_var]]
    pooled_forecast_errors[[fold]] = data.frame(fe = errors, class = classes) %>%
      tibble::rownames_to_column(var = "obs_no")
    
  }
  
  # looks like our predictions are aligned with our observations
  
  # forecast_errors[[5]] %>%
  #   ggplot2::ggplot(ggplot2::aes(x = y, y = yhat)) + 
  #   ggplot2::geom_point() + 
  #   ggplot2::facet_grid(class ~ knot) + 
  #   ggplot2::geom_abline()
  
  # OK now we have the forecast_error list for each fold.  Let's combine it 
  # into a single data frame
  
  fe = dplyr::bind_rows(forecast_errors, .id = "fold") %>% 
    dplyr::mutate(knot = as.numeric(knot))
  
  pooled_fe = dplyr::bind_rows(pooled_forecast_errors, .id = "fold") 
  
  pooled_class_fold_summary = pooled_fe %>% 
    dplyr::group_by(class, fold) %>% 
    dplyr::summarise(
      MSFE = mean(fe^2, na.rm = TRUE),
      MAFE = mean(abs(fe), na.rm = TRUE)      
    ) %>% 
    dplyr::ungroup()
  
  pooled_class_summary = pooled_class_fold_summary %>% 
    dplyr::group_by(class) %>% 
    dplyr::summarise(
      MSFE = mean(MSFE),
      MAFE = mean(MAFE)
    ) %>% 
    dplyr::ungroup()
  
  
  pooled_summary = pooled_class_summary %>% 
    dplyr::summarise(
      MSFE = mean(MSFE),
      MAFE = mean(MAFE)
    )
  
  # now we need to work out the performance measures for each knot
  # within each class, knot and fold, calculate performance
  class_knot_fold_summary = fe %>% 
    dplyr::group_by(class, knot, fold) %>% 
    dplyr::summarise(
      MSFE = mean(fe^2, na.rm = TRUE),
      MAFE = mean(abs(fe), na.rm = TRUE)
    ) %>% 
    dplyr::ungroup()
  
  # class_knot summary (i.e. averaging over the K folds)
  class_knot_summary = class_knot_fold_summary %>% 
    dplyr::group_by(class, knot) %>% 
    dplyr::summarise(
      MSFE = mean(MSFE),
      MAFE = mean(MAFE)
    ) %>% 
    dplyr::ungroup()
  
  
  # knot summary (ignores class size)
  knot_summary = class_knot_summary %>% 
    dplyr::group_by(knot) %>% 
    dplyr::summarise(
      MSFE = mean(MSFE),
      MAFE = mean(MAFE)
    ) %>% 
    dplyr::ungroup()
  
  
  # "best" knot
  
  knot_MSFE = knot_summary %>% 
    dplyr::arrange(MSFE) %>% 
    dplyr::slice(1) %>%
    dplyr::pull(knot)
  knot_MAFE = knot_summary %>% 
    dplyr::arrange(MAFE) %>% 
    dplyr::slice(1) %>%
    dplyr::pull(knot)
  
  # identify and return "best" fits
  # missing last knot?!?
  fit_MSFE = fit$coef_df %>% dplyr::filter(knot == knot_MSFE)
  fit_MAFE = fit$coef_df %>% dplyr::filter(knot == knot_MAFE)
  

  
  return(list(train = train, test = test, cv_fits = cv_fits,
              fe = fe, pooled_fe = pooled_fe,
              class_knot_fold_summary = class_knot_fold_summary,
              class_knot_summary = class_knot_summary, 
              knot_summary = knot_summary,
              pooled_class_fold_summary = pooled_class_fold_summary,
              pooled_class_summary = pooled_class_summary,
              pooled_summary = pooled_summary,
              fit_MSFE = fit_MSFE,
              fit_MAFE = fit_MAFE))
  
}


#' Function to compute optimal value of regularization parameter based on K-fold cross-validation
#' @param fit fitted object returned from multiclass_reg()
#' @param Kfold K-fold cross-validation
#' @param nweight vector of length K indicating weights for MSFE measure
#' @param asym flag: asymmetric forecast error measure or not
#' @param type MSFE or MAFE
#' @param l1length length sparsity grid l1 penalty
#' @return A list with the following components
#' \item{lambda_opt}{Selected value of the regularization parameter via BIC}
#' \item{bhat_opt}{Estimated coefficients selected via BIC}
#' @examples
#' p <- 7
#' k <- 2
#' n <- 50
#' beta <- c(1, 2)
#' set.seed(1)
#' X <- list(matrix(rnorm(p * n), ncol = p), matrix(rnorm(p * n), ncol = p))
#' Y <- list(rnorm(n), rnorm(n))
#' res <- multiclass_reg(Y, X)
#' cv_res <- multiclass_cv(res)
#' @export
multiclass_cv <- function(fit, # bhats, lambdas,  Y, X, extract these from the fitted object
                          Kfold = 5,
                          nweight = NULL,
                          weighting = TRUE, wtype = "size",
                          type = "MSFE", l1length = 100,
                          normalize = TRUE) {
  
  bhats <- fit$coef
  lambdas <- fit$lambdas
  Y <- fit$Y_list
  X <- fit$X_list
  
  p <- unique(sapply(X, ncol))
  K <- length(Y) # number of classes
  
  determine.interval <- function(FE, lfknot) {
    # MSFEcv needs to be of dimension Kfold x nbrknots
    MSFEavg <- apply(FE, 2, mean, na.rm = T)
    MSFEsd <- apply(FE, 2, sd, na.rm = T)
    ub <- MSFEavg[lfknot] + MSFEsd[lfknot] / sqrt(dim(FE)[1])
    MSFEup <- apply(FE, 2, mean, na.rm = T)
    MSFEup[1:lfknot] <- NA
    knotub <- which.min(abs(MSFEup - ub))
    MSFElb <- apply(FE, 2, mean, na.rm = T)
    MSFElb[-c(1:lfknot)] <- NA
    knotlb <- which.min(abs(MSFElb - ub))
    
    out <- list("knotUB" = knotub, "knotLB" = knotlb)
  }
  
  determine.weights <- function(FE, knots) {
    # FE needs to be a vector of length nbr.knots
    
    # MSFE-based weights
    wsel <- (FE[min(knots):max(knots)]^-1) / sum((FE[min(knots):max(knots)]^-1))
    knotp <- min(which(FE == FE[length(FE)]))
    wall <- (FE[1:knotp]^-1) / sum((FE[1:knotp]^-1))
    
    # rank-based weights
    rwsel <- (rank(FE[min(knots):max(knots)])^-1) / sum((rank(FE[min(knots):max(knots)])^-1))
    rwall <- (rank(FE[1:knotp])^-1) / sum((rank(FE[1:knotp])^-1))
    
    out <- list("wsel" = wsel, "wall" = wall, "rwsel" = rwsel, "rwall" = rwall, "knotp" = knotp)
  }
  
  # Store Forecast Error Results
  MSFEcv <- array(NA, c(Kfold, length(lambdas), l1length))
  MSFEKcv <- array(NA, c(length(lambdas), K, Kfold, l1length))
  
  # Preliminaries : set dimensions
  K <- length(Y)
  nseq <- sapply(Y, length)
  if (is.null(nweight)) {
    # nweight <- rank(nseq)^(-1)/sum(rank(nseq)^(-1))
    nweight <- rep(1 / K, K)
  }
  
  # Observations in random order
  nrandom <- vector("list", K)
  for (i.K in 1:K) {
    # nrandom[[i.K]] <- sample(1:nseq[i.K], nseq[i.K])
    nrandom[[i.K]] <- 1:nseq[i.K]
  }
  
  # Start cross-validation loop
  for (i.cv in 1:Kfold) {
    # cat("start", i.cv, "\n")
    Ytrain <- Ytest <- Xtrain <- Xtest <- vector("list", K)
    ind_group <- c()
    
    for (i.K in 1:K) {
      cv_cut <- ceiling(nseq[i.K] / Kfold)
      ind_group <- c(ind_group, rep(i.K, cv_cut))
      test_set <- nrandom[[i.K]][(1 + (cv_cut) * (i.cv - 1)):(cv_cut * i.cv)]
      
      
      if (all(is.na(test_set))) {
        Ytrain[[i.K]] <- Y[[i.K]]
        Xtrain[[i.K]] <- X[[i.K]]
        
        Ytest[[i.K]] <- rep(NA, cv_cut)
        Xtest[[i.K]] <- matrix(999, nrow = cv_cut, ncol = p)
        Xtest[[i.K]][, is.na(colSums(Xtrain[[i.K]]))] <- NA
      } else {
        if (any(is.na(test_set))) {
          mark.na <- which(is.na(test_set))
          Ytrain[[i.K]] <- Y[[i.K]][-test_set[-mark.na]]
          Xtrain[[i.K]] <- X[[i.K]][-test_set[-mark.na], ]
          
          Ytest[[i.K]] <- c(Y[[i.K]][test_set[-mark.na]], rep(NA, length(mark.na)))
          Xtest[[i.K]] <- matrix(999, nrow = cv_cut, ncol = p)
          Xtest[[i.K]][, is.na(colSums(Xtrain[[i.K]]))] <- NA
          Xtest[[i.K]][1:(cv_cut - length(mark.na)), ] <- X[[i.K]][test_set[-mark.na], ]
        } else {
          Ytrain[[i.K]] <- Y[[i.K]][-test_set]
          Xtrain[[i.K]] <- X[[i.K]][-test_set, ]
          
          Ytest[[i.K]] <- Y[[i.K]][test_set]
          Xtest[[i.K]] <- X[[i.K]][test_set, ]
          
          if (cv_cut == 1) {
            Xtest[[i.K]] <- matrix(Xtest[[i.K]], nrow = 1)
          }
        }
      }
      colnames(Xtest[[i.K]]) <- colnames(X[[1]])
    }
    
    mreg <- multiclass_reg(
      Y = Ytrain,
      X = Xtrain,
      weighting = weighting,
      wtype = wtype,
      normalize = normalize
    )
    
    # coefs <- matrix(NA, ncol=length(lambdas), nrow=p*K) #  USING p NOT OK!!!
    coefs <- matrix(NA, ncol = length(lambdas), nrow = dim(mreg$coef)[1]) #
    if (min(lambdas) < min(mreg$fit$lambda) && !mreg$fit$completepath) {
      index <- which(lambdas < min(mreg$fit$lambda))
      beta <- c()
      for (it in 1:K) {
        Xk <- Xtrain[[it]]
        Xk <- Xk[, !is.na(colSums(Xk))]
        Yk <- matrix(Ytrain[[it]], ncol = 1)
        beta <- c(beta, solve(t(Xk) %*% Xk + 0.001) %*% t(Xk) %*% Yk) # Small ridge penalty if LS not computable
      }
      coefs[, index] <- beta
      coefs[, -index] <- coef(mreg$fit, lambdas[-index])$beta
    } else {
      coefs <- coef(mreg$fit, lambdas)$beta
    }
    
    # Sparse Solution
    l1pen <- sort(c(0, exp(seq(log(max(abs(coefs))), log(min(abs(coefs))), length = l1length - 1))))
    
    shrinkcoefs <- array(NA, c(dim(coefs), length(l1pen)))
    
    testfit <- inputs(
      Y = Ytest,
      X = Xtest,
      buildD = FALSE
    )
    
    if (length(which(testfit$X == 999)) != 0) {
      testfit$X[which(testfit$X == 999)] <- NA
    }
    
    Yhat <- vector("list", length(l1pen))
    
    for (il1 in 1:l1length) {
      shrinkcoefs[, , il1] <- apply(coefs, 2,
                                    lambda = l1pen[il1],
                                    function(U, lambda) {
                                      sign(U) * pmax(abs(U) - lambda, 0)
                                    }
      )
      Yhat[[il1]] <- testfit$X %*% shrinkcoefs[, , il1]
    }
    
    # Forecast Errors for each group separately
    MSFEgroup <- array(NA, c(ncol(Yhat[[1]]), K, length(l1pen)))
    
    for (i.k in 1:K) {
      groupiK <- which(ind_group == i.k)
      nbr.test <- sapply(Ytest, length)[i.k]
      
      if (nbr.test == 1) {
        if (type == "MSFE") {
          for (i.l1 in 1:l1length) {
            MSFEgroup[, i.k, i.l1] <- (testfit$Y[groupiK] - Yhat[[i.l1]][groupiK, ])^2
          }
        }
        
        if (type == "MAFE") {
          for (i.l1 in 1:l1length) {
            MSFEgroup[, i.k, i.l1] <- abs(testfit$Y[groupiK] - Yhat[[i.l1]][groupiK, ])
          }
        }
      } else {
        if (type == "MSFE") {
          for (i.l1 in 1:l1length) {
            MSFEgroup[, i.k, i.l1] <- apply((testfit$Y[groupiK] - Yhat[[i.l1]][groupiK, ])^2,
                                            2, mean,
                                            na.rm = T
            )
          }
        }
        
        if (type == "MAFE") {
          for (i.l1 in 1:l1length) {
            MSFEgroup[, i.k] <- apply(abs(testfit$Y[groupiK] - Yhat[[i.l1]][groupiK, ]),
                                      2, mean,
                                      na.rm = T
            )
          }
        }
      }
    }
    
    MSFEKcv[, , i.cv, ] <- MSFEgroup
    
    # Forecast Errors averaged over the different groups
    for (i.l1 in 1:l1length) {
      MSFEcv[i.cv, , i.l1] <- apply(MSFEgroup[, , i.l1], 1, Kweight = nweight, function(U, Kweight) {
        sum(U * Kweight, na.rm = T)
      })
    }
  }
  
  #### Overall MSFE curve, averaged over all knots ####
  MSFEavg <- apply(MSFEcv, c(2, 3), mean, na.rm = T) # Of Dimension (knots)x(l1length)
  # Choose best value for l1pen
  l1opt.knot <- which.min(apply(MSFEavg, 2, min))
  
  # Optimal "fusion knot" at optimal l1-MSFE curve
  lfl1opt.knot <- which.min(MSFEavg[, l1opt.knot])
  lfl0.knot <- which.min(MSFEavg[, 1])
  
  # Determine interval around min knot
  lfl1opt.lbub <- unlist(determine.interval(FE = MSFEcv[, , l1opt.knot], lfknot = lfl1opt.knot))
  lfl0.lbub <- unlist(determine.interval(FE = MSFEcv[, , 1], lfknot = lfl0.knot))
  
  # Determine Forecast Combination Weights
  lfl0.FC <- determine.weights(FE = MSFEavg[, 1], knots = lfl0.lbub)
  lfl1opt.FC <- determine.weights(FE = MSFEavg[, l1opt.knot], knots = lfl1opt.lbub)
  
  #### MSFE curve per group ####
  MSFEKavg <- apply(MSFEKcv, c(1, 2, 4), mean, na.rm = T) # of dimension nbr.knots x K x l1length
  l1opt.knotK <- lfl1opt.knotK <- lfl0.knotK <- rep(NA, K)
  lfl1opt.lbubk <- lfl0.lbubk <- matrix(NA, ncol = 2, nrow = K)
  colnames(lfl1opt.lbubk) <- colnames(lfl0.lbubk) <- c("LBk", "UBk")
  lfl0.FCk <- lfl1opt.FCk <- vector("list", K)
  for (ik in 1:K) {
    l1opt.knotK[ik] <- which.min(apply(MSFEKavg[, ik, ], 2, min)) # Optimal l1 knot for each class
    lfl1opt.knotK[ik] <- which.min(MSFEKavg[, ik, l1opt.knotK[ik]]) # Optimal "fusion" knot at optimal l1 curve
    lfl0.knotK[ik] <- which.min(MSFEKavg[, ik, 1]) # Optimal "fusion" knot at l1=0 curve
    
    lfl0.lbubk[ik, ] <- unlist(determine.interval(FE = t(MSFEKcv[, ik, , 1]), lfknot = lfl0.knotK[ik]))
    lfl1opt.lbubk[ik, ] <- unlist(determine.interval(FE = t(MSFEKcv[, ik, , l1opt.knotK[ik]]), lfknot = lfl1opt.knotK[ik]))
    
    lfl0.FCk[[ik]] <- determine.weights(FE = apply(MSFEKcv[, ik, , 1], 1, mean, na.rm = T), knots = lfl0.lbubk[ik, ])
    lfl1opt.FCk[[ik]] <- determine.weights(FE = apply(MSFEKcv[, ik, , l1opt.knotK[ik]], 1, mean, na.rm = T), knots = lfl1opt.lbubk[ik, ])
  }
  
  
  # MSFEcv : of dimension Kfold x nbr.knots x l1length
  # MSFEavg : of dimension nbr.knots x l1length
  # l1opt.knot : optimal value sparsity parameter (MSFE averaged over all groups)
  # lfl1opt.knot : optimal fusion parameter at optimal l1 curve
  # lfl0.knot : optimal fusion parameter at l1=0 curve
  # lfl1opt.lbub : LB and UB interval around lfl1opt
  # lfl0.lbub : LB and UB interval around lfl0
  # lfl0.FC : list with FC results for lfl0
  # lfl1opt.FC : list with FC results for lfl1opt
  # MSFEKcv : of dimension nbr.knots x K x Kfold x l1length
  # MSFEKavg : of dimension nbr.knots x K x l1length
  # l1opt.knotK : vector optimal value sparsity parameter for each group
  # lfl1opt.knotK : vector optimal value fusion parameter at optimal l1 curve for each group
  # lfl0.knotK : vector optimal value sparsity parameter at l1=0 curve each group
  # lfl0.lbubk : K x 2 matrix : LB and UB around lfl0k
  # lfl1opt.lbubk : K x 2 matrix : LB and UB around lfl1optk
  # lfl0.FCk : list with FC results for each knot k for lfl0
  # lfl1opt.FCk : list with FC results for each knot k for lfl1opt
  
  out <- list(
    "MSFEcv" = MSFEcv, 
    "MSFEavg" = MSFEavg,
    "l1opt.knot" = l1opt.knot, 
    "lfl1opt.knot" = lfl1opt.knot, 
    "lfl0.knot" = lfl0.knot,
    "lfl1opt.lbub" = lfl1opt.lbub, 
    "lfl0.lbub" = lfl0.lbub, 
    "lfl0.FC" = lfl0.FC, 
    "lfl1opt.FC" = lfl1opt.FC,
    "MSFEKcv" = MSFEKcv, 
    "MSFEKavg" = MSFEKavg,
    "l1opt.knotK" = l1opt.knotK, 
    "lfl1opt.knotK" = lfl1opt.knotK, "lfl0.knotK" = lfl0.knotK,
    "lfl0.lbubk" = lfl0.lbubk, "lfl1opt.lbubk" = lfl1opt.lbubk, "lfl0.FCk" = lfl0.FCk, "lfl1opt.FCk" = lfl1opt.FCk
  )
}
