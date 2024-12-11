#' svThresholdOp
#'
#' @description
#' Implements the Singular Value Thresholding operator as described in 'A SINGULAR VALUE THRESHOLDING ALGORITHM FORMATRIX COMPLETION', Cai et al. (2010)
#'
#'
#' @param M A matrix
#' @param tau A scalar, used for soft shrinkage of singular values
#'
#' @return Matrix with same dimensions as M, shrunk according to the operator.
svThresholdOp <- function(M, tau){
  #Get svd
  svM <- svd(M)

  #Shrink SV's
  svShrink <- pmax(svM$d - tau, rep(0, length(svM$d)))

  return(svM$u %*% diag(svShrink) %*% t(svM$v))
}
