#' Eigenvector Clustering with Optional Thresholding and Minimum Cluster Size
#'
#' This function performs eigenvector clustering on a given symmetric (co-association) matrix.
#' The function can optionally normalize the input matrix, and use the threshold method (if given) to
#' break up clusters that should not be attached. The minimum allowable cluster size can
#' also be specified.
#'
#' @param A A symmetric \code{n x n} matrix on which clustering is
#' to be performed. Assumed to be a type of co-assocation matrix. That is,
#' higher values in \eqn{A_{ij}} imply that observations \eqn{i} and \eqn{j} are **more** likely
#' to be in the same cluster.
#' @param threshold A numeric value between 0 and 1 (or NULL). This parameter controls
#' the average strength of association points must have with the eigenvector in order to NOT break the
#' cluster down further. If the average strength is lower than the given threshold the cluster,
#' will be broken down by next largest eigenvector and sign. If NULL, this step is skipped.
#' Defaults to NULL.
#' @param minSize An integer specifying the minimum allowed cluster size. Defaults to 1.
#' @param normalize A logical value indicating whether to normalize \eqn{\Big(A_{Norm} = D_A^{-1/2}AD_A^{-1/2}\Big)} the matrix
#' before clustering. Defaults to \code{FALSE}.
#'
#' @return A numeric vector of length \code{n} representing the cluster
#' assignment of each point.
#'
#' @details
#' The function first checks that the input matrix \code{A} is square and symmetric.
#' If the \code{normalize} argument is set to \code{TRUE}, the function normalizes
#' the matrix by scaling it using the diagonal values. If \code{replicate} is not
#' zero, a portion of the matrix is replicated. The function then performs an
#' eigendecomposition of the matrix and uses the resulting eigenvectors and
#' eigenvalues to cluster the data. If the \code{minSize} argument is greater than 1,
#' the function ensures that no cluster is smaller than the specified minimum size
#' by reassigning points from smaller clusters to the next closest option based on the eigendecomposition.
#'
#' @note The \code{replicate} value should typically be between 0 and 1.
#' Setting \code{replicate} to a value greater than or equal to 1 may lead
#' to unexpected behavior.
#'
#' @examples
#' # Create a symmetric matrix
#' A <- matrix(runif(25), nrow = 5)
#' A <- (A + t(A)) / 2
#'
#' # Perform clustering with default parameters and no thresholding
#' clust <- eCl.threshold(A)
#'
#' # Perform clustering with suggested threshold and minimum cluster size
#' clust <- eCl.threshold(A, threshold = sqrt(.5), minSize = 2)
#'
#' @export
eCl.threshold <- function(A,
                          threshold = NULL,
                          minSize = 1,
                          normalize = FALSE){
  #Test to make sure A is as expected
  if(!is.matrix(A)){
    stop('A must be a n x n matrix.')
  }
  if(nrow(A) != ncol(A)){
    stop('A must be a n x n matrix.')
  }
  if(!isSymmetric(A)){
    stop('A must be symmetric')
  }

  #Test to make sure threshold is reasonable
  if(!is.null(threshold)){
    if(threshold < 0){
      stop('threshold cannot be negative')
    }
    if(threshold>= 1){
      warning('threshold value of >= 1 is not allowed')
    }
  }

  #Get sample size
  n <- nrow(A)

  #Make sure minimum cluster size makes any sense
  if(minSize > n){
    stop('minSize cannot be larger than the sample size (number of rows of A)')
  }

  #Normalize if requested
  if(normalize){
    dA <- diag(sqrt(1/diag(A)))
    A <- dA %*% A %*% dA
  }

  #Eigendecomposition
  eDec <- eigen(A, symmetric = TRUE)
  eVec <- eDec$vectors
  eVal <- eDec$values
  eVal <- pmax(eVal, rep(0, length(eVal)))

  #Get L
  L <- eVec %*% diag(sqrt(eVal))

  #Get clustering
  clust <- sapply(1:n, function(i){
    Li <- abs(L[i,])
    return(which(Li == max(Li))[1])
  })

  #If minSize > 1, check and redo if necessary
  if(minSize > 1){
    cSizes <- table(clust)
    t <- rep(2,n) #Initialize which loop counter
    minClust <- min(cSizes) #initialize minCluster size
    while(minClust < minSize){
      #Identify smallest cluster(s)
      tooSmall <- names(cSizes)[as.numeric(which(cSizes == minClust))]

      #Identify elements in tooSmall clusters
      tooSmalli <- which(clust %in% tooSmall)

      #Update for each element
      for(i in tooSmalli){
        ti <- t[i]
        Li <- abs(L[i,])
        sortedLi <- sort(Li, decreasing = TRUE)
        clust[i] <- which(Li == sortedLi[ti])[1]
        t[i] <- ti + 1
      }

      #Update loop parameters
      cSizes <- table(clust)
      minClust <- min(cSizes)
    }
  }

  #If threshold is NOT null, perform thresholding procedure
  if(!is.null(threshold)){
    #all clusters
    allClust <- unique(clust)

    #Now, for each cluster with average strength under threshold
    #split based on next SIGN of next largest eigenvector
    subClusts <- lapply(allClust, function(cl){
      #Get cluster specific A
      A.cl <- A[clust == cl, clust == cl]

      #Eigendecomposition for cluster
      eDec.cl <- eigen(A.cl, symmetric = TRUE)
      eVec.cl <- eDec.cl$vectors
      eVal.cl <- eDec.cl$values
      eVal.cl <- pmax(eVal.cl, rep(0, length(eVal.cl)))

      #Get new L
      L.cl <- eVec.cl %*% diag(sqrt(eVal.cl))

      lCl <- L.cl[,1]

      #Check each cluster for average strength with eigenvector
      avgStrength <- mean(abs(lCl))

      if(avgStrength < threshold){
        #Pick the column from L with the SECOND largest |colMean|
        #L2 <- L.cl[,which(colMeans(abs(L.cl)) == sort(colMeans(abs(L.cl)),
                                                      #decreasing = TRUE)[2])]
        L2 <- L.cl[,2]

        #Sort by sign
        newClust <- (sign(L2) + 3)/2 #+3/2 to make it 1/2 instead of -1/1

        #Check to make sure it passes minSize requirements
        if(all(table(newClust) >= minSize)){
          return(newClust)
        }else{
          return(rep(1, sum(clust == cl)))
        }
      }else{
        return(rep(1, sum(clust == cl)))
      }
    })

    clust2 <- clust
    for(i in 1:length(allClust)){
      clust2[clust == allClust[i]] <- subClusts[[i]]
    }
    clust <- interaction(clust, clust2)
  }

  return(clust)
}
