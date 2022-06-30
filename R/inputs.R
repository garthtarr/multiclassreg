#' Function to get inputs for fused lasso fit
#'
#' @param Y list of size K. Its kth entry is a vector of size n_k
#' @param X list of size K. Its kth entry is a matrix of size n_kxp
#' @param weighting logical. If TRUE: use weighted regularization parameter
#' @param weights provide vector of weights
#' @param wtype either `"size"` (default, weighting based on group sample sizes) or `"sizesqrt"`
#' @param buildD logical. If TRUE: compute Dmatrix
#' @param normalize logical. Whether or not to normalize the weights (currently divides by maximum).
#' @return A list with the following components
#' \item{Y}{Response vector}
#' \item{X}{Block diagonal design matrix}
#' \item{D}{Coupling matrix}
#' \item{var_indicator}{pxK : indicates corresponding column (in the D-matrix) of each variable in each class}
#' @export
inputs <- function(
  Y, X, 
  weighting = FALSE, 
  weights = NULL,
  wtype = "size",
  buildD = TRUE, 
  normalize = TRUE
  ){
  
  # Start code
  # Dimensions
  K <- length(Y) # number f classes
  # Effective number of predictors per class
  pseq <- sapply(X, function(U){ sum(!is.na(apply(U, 2, mean, na.rm=T))) }) 
  p = ncol(X[[1]])
  n <- sapply(Y, length)
  
  
  var_indicator <- NA
  Dmatrix <- NA
  if(buildD){
    
    # Prelim to buildD D-matrix when the number of variables varies across classes
    # var_indicator contains in cell i,j the column number in the D-matrix of variable i in class j
    var_indicator <- matrix(NA, ncol = K, nrow = p)
    count <- 0
    
    for(i.K in 1:K){
      for(i.p in 1:p){ 
        if(!any(is.na(X[[i.K]][, i.p]))){ 
          count <- count + 1
          var_indicator[i.p, i.K] <- count
        }
      }
    }
    
    # D-matrix that represents the coupling
    Dmatrix <- c()
    for(i.p in 1:p){
      for(i.K1 in 1:(K-1)){
        for(i.K2 in (i.K1 +1) :K){
          if(!is.na(var_indicator[i.p, i.K1]) & !is.na(var_indicator[i.p, i.K2])){
            new <- matrix(0, ncol = max(var_indicator, na.rm = T), nrow=1)
            new[1, c(var_indicator[i.p, i.K1], var_indicator[i.p, i.K2])] <- c(1, -1)
            Dmatrix <- rbind(Dmatrix, new)
          }
        }
      }
    }
    
    if(weighting){
      
      if(is.null(weights)){ # No weights provided
        weights <- rep(NA, nrow(Dmatrix)) # weighting vector
      }
      
      if(wtype =="size" | wtype == "sizesqrt"){
        ratios <- c()
        for(i.D in 1:nrow(Dmatrix)){
          index <- 1:K
          classseq <- unlist(lapply(index, function(u, pseq){rep(u, pseq[u])}, pseq=pseq)) # indicates to which class each column of Dmatrix belongs
          sizes <- n[classseq[which(Dmatrix[i.D,]!=0)]]
          
          if(wtype=="size"){
            ratios <- c(ratios, max(sizes)/min(sizes))
            weights[i.D] <- max(sizes)/min(sizes)
          }
          
          if(wtype=="sizesqrt"){
            ratios <- c(ratios, max(sizes)/min(sizes))
            weights[i.D] <- sqrt(max(sizes)/min(sizes))
          }
        }
      }
      # print(weights)
      if(normalize==T){
        weights <- weights/max(ratios)
      }
      
      # print(weights)
      Dmatrix <- Dmatrix*weights
    }
    
  }
  
  # Build predictor set deleting the NA columns
  Xnew <- vector("list", K)
  for(i.K in 1:K){
    Xnew[[i.K]] <- X[[i.K]][, which(!is.na(var_indicator[, i.K]))]
    colnames(Xnew[[i.K]]) <- colnames(X[[i.K]])[which(!is.na(var_indicator[, i.K]))]
  }
  
  # Create Y and X inputs for fusedlasso function
  Ypool <- unlist(Y)
  Xkron <- as.matrix(Matrix::bdiag(Xnew))
  colnames(Xkron) = unlist(lapply(Xnew, colnames))
  
  out <- list("Y" = Ypool, "X" = Xkron, "D" = Dmatrix, "var_indicator" = var_indicator)
}

