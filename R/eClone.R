#' Eigenvector Clustering with Optional Replication and Minimum Cluster Size
#'
#' This function performs eigenvector clustering on a given symmetric (co-association) matrix.
#' The function can optionally normalize the input matrix, and replicate a
#' a single value in the matrix before clustering. The minimum allowable cluster size can
#' also be specified.
#'
#' @param A A symmetric \code{n x n} matrix on which clustering is
#' to be performed. Assumed to be a type of co-assocation matrix. That is,
#' higher values in \eqn{A_{ij}} imply that observations \eqn{i} and \eqn{j} are **more** likely
#' to be in the same cluster.
#' @param replicate A numeric value between 0 and 1. This parameter controls
#' the extent of replication of a selected point in the matrix. If nonzero a random
#' point is replicated \code{ceiling(replicate*n)} times. Defaults to 0.
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
#' # Perform clustering with default parameters
#' clust <- eCl.one(A)
#'
#' # Perform clustering with replication and minimum cluster size
#' clust <- eCl.one(A, replicate = 0.1, minSize = 2)
#'
#' @export
eCl.one <- function(A,
                    replicate = 0,
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

  #Test to make sure replicate is reasonable
  if(replicate < 0){
    stop('replicate cannot be negative')
  }
  if(replicate >= 1){
    warning('replication value of >= 1 is probably unreasonable and will break things')
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

  #If replication, replicate
  if(replicate != 0){
    #pick a point to replicate
    a <- sample(n, 1)

    #replication number is ceiling of sampleSize*replicate
    repNum <- ceiling(n*replicate)

    #1 grid
    g1 <- matrix(1, nrow = repNum, ncol = repNum)

    #bottom
    botA <- A[rep(a,repNum),, drop = F]

    #side
    sideA <- A[,rep(a,repNum), drop = F]

    #Attach
    A <- cbind(rbind(A, botA), rbind(sideA,g1))
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

  return(clust)
}
