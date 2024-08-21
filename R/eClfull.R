#' Full Eigenvector Clustering with Optional Replication and Minimum Cluster Size
#'
#' This function performs a full eigenvector clustering on a given symmetric (co-association) matrix
#' with an initial pass followed by a phase focused on separating the initial clusters if needed. It allows for an
#' optional minimum cluster size. This function builds upon the `eCl.one` function.
#'
#' @inheritParams eCl.one
#'
#' @return A factor vector of length \code{n} representing the final cluster
#' assignment of each point after the full clustering process.
#'
#' @details
#' The function first normalizes the input matrix \code{A} if requested.
#' It then performs an initial pass of clustering using the \code{\link{eCl.one}}
#' function with replication turned off. During the breaking phase, the function
#' refines the clusters by applying \code{\link{eCl.one}} to each unique cluster
#' with the specified \code{replicate} value. The final cluster assignments are
#' obtained by combining the results of the initial and breaking phases.
#'
#' @note The \code{replicate} value should typically be between 0 and 1.
#' Setting \code{replicate} to a value greater than or equal to 1 may lead
#' to unexpected behavior.
#'
#' @seealso \code{\link{eCl.one}} for the underlying clustering function.
#'
#' @examples
#' # Create a symmetric matrix
#' A <- matrix(runif(25), nrow = 5)
#' A <- (A + t(A)) / 2
#'
#' # Perform full clustering with default parameters
#' clust <- eCl.full(A)
#'
#' # Perform full clustering with normalization and minimum cluster size
#' clust <- eCl.full(A, normalize = TRUE, minSize = 2)
#'
#' @export
eCl.full <- function(A,
                     replicate = .2,
                     minSize = 1,
                     normalize = FALSE){
  #Normalize if requested
  if(normalize){
    dA <- diag(sqrt(1/diag(A)))
    A <- dA %*% A %*% dA
  }

  #Run first pass (no replication allowed on first pass)
  clust <- eCl.one(A = A, replicate = 0, minSize = minSize)

  ########################################
  #Break Phase
  ########################################
  uniqueC <- unique(clust) #Unique clusters

  #Break per unique cluster
  subClusts <- lapply(uniqueC, function(uC){
    subA <- A[clust == uC, clust == uC, drop = F]

    return(eCl.one(A = subA,
                   replicate = replicate,
                   minSize = minSize))
  })

  #Build interaction of clusters
  clust2 <- clust
  for(i in 1:length(uniqueC)){
    clust2[clust == uniqueC[i]] <- subClusts[[i]]
  }

  clust <- interaction(clust, clust2)

  return(clust)
}
