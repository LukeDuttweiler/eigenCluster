#' Eigenvector Clustering with Optional Minimum Cluster Size
#'
#' This function performs eigenvector clustering on a given symmetric (co-association) matrix.
#' The function can optionally normalize the input matrix, and uses the reflection method to
#' break up clusters that should not be attached. The minimum allowable cluster size can
#' also be specified.
#'
#' @param A A symmetric \code{n x n} matrix on which clustering is
#' to be performed. Assumed to be a type of co-assocation matrix. That is,
#' higher values in \eqn{A_{ij}} imply that observations \eqn{i} and \eqn{j} are **more** likely
#' to be in the same cluster.
#' @param minSize An integer specifying the minimum allowed cluster size. Defaults to 1.
#' @param normalize A logical value indicating whether to normalize \eqn{\Big(A_{Norm} = D_A^{-1/2}AD_A^{-1/2}\Big)} the matrix
#' before clustering. Defaults to \code{FALSE}.
#' @param maxK An integer \code{<= n} specifying that maximum possible number of clusters in the data.
#' Speeds computation. Default is min(100, n).
#'
#' @return A numeric vector of length \code{n} representing the cluster
#' assignment of each point.
#'
#' @details
#' The function first checks that the input matrix \code{A} is square and symmetric.
#' If the \code{normalize} argument is set to \code{TRUE}, the function normalizes
#' the matrix by scaling it using the diagonal values. The function then performs an
#' eigendecomposition of the matrix and uses the resulting eigenvectors and
#' eigenvalues to cluster the data, breaking and building clusters as specified
#'  by the break-build eigencluster method. If the \code{minSize} argument is greater than 1,
#' the function ensures that no cluster is smaller than the specified minimum size
#' by reassigning points from smaller clusters to the next closest option based on the eigendecomposition.
#' maxK allows the user to set an argument for maximum possible clusters. For datasets with a large
#' sample size this allows a much faster computation time. The default maxK is 100, and probably doesn't
#' need to be set much smaller than that in general.
#'
#'
#' @examples
#' # Create a symmetric matrix
#' A <- matrix(runif(25), nrow = 5)
#' A <- (A + t(A)) / 2
#'
#' # Perform clustering with default parameters and no thresholding
#' clust <- eCl.breakBuild(A)
#'
#' # Perform clustering with minimum cluster size
#' clust <- eCl.breakBuild(A, minSize = 2)
#'
#' @export
eCl.breakBuild <- function(A,
                           minSize = 1,
                           normalize = FALSE,
                           maxK = min(100,nrow(A))){
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

  ###############
  #INITIAL BREAK
  ###############

  if(nrow(A) > maxK){
    eDec <- Rfast::eigen.sym(A, k = min(maxK, nrow(A)))
    eVec <- eDec$vectors
    eVal <- eDec$values[,1]
  }else{
    eDec <- eigen(A, symmetric = TRUE)
    eVec <- eDec$vectors
    eVal <- eDec$values
  }

  eVal <- pmax(eVal, rep(0, length(eVal)))

  #Get L
  L <- eVec %*% diag(sqrt(eVal))

  #Get clustering
  clust <- apply(abs(L), 1, which.max)
  clust <- as.character(clust)

  ############
  #SUB-BREAK 1
  ############
  allClust <- unique(clust)
  done <- rep(FALSE, length(allClust))
  t <- 1

  #Now, for each cluster split based on eigendecomp of principal submatrix until no longer works
  while(!all(done)){
    cl <- allClust[t]
    #Get cluster specific A
    A.cl <- A[clust == cl, clust == cl, drop = F]

    #As long as cluster isn't singleton, try to break down
    if(nrow(A.cl) > 1){
      if(nrow(A.cl) > maxK){
        eDec.cl <- Rfast::eigen.sym(A.cl, k = min(maxK, nrow(A.cl)))
        eVec.cl <- eDec.cl$vectors
        eVal.cl <- eDec.cl$values[,1]
      }else{
        eDec.cl <- eigen(A.cl, symmetric = TRUE)
        eVec.cl <- eDec.cl$vectors
        eVal.cl <- eDec.cl$values
      }

      eVal.cl <- pmax(eVal.cl, rep(0, length(eVal.cl)))
      L.cl <- eVec.cl %*% diag(sqrt(eVal.cl))

      whichL <- apply(abs(L.cl), 1, which.max) #Which L is each obs in cluster cl assocaiated with?
      newClust <- whichL
    }else{
      newClust <- rep(1, nrow(A.cl))
    }

    clust[clust == cl] <- as.character(interaction(cl, newClust))

    #Check to see if done
    if(!all(newClust == 1)){
      allClust <- c(allClust, unique(as.character(interaction(cl, newClust))))
      done <- c(done, rep(FALSE, length(unique(newClust))))
    }

    done[t] <- TRUE
    t <- t+1
  }

  #subClusts <- lapply(allClust, function(cl){

    #Get cluster specific A
  #  A.cl <- A[clust == cl, clust == cl, drop = F]

    #As long as cluster isn't singleton, try to break down
  #  if(nrow(A.cl) > 2){
      #Eigendecomposition for cluster
      #eDec.cl <- suppressWarnings(RSpectra::eigs_sym(A.cl, k = min(maxK, nrow(A.cl))))
      #eVec.cl <- eDec.cl$vectors
      #eVal.cl <- eDec.cl$values
      #eVal.cl <- pmax(eVal.cl, rep(0, length(eVal.cl)))

      #Get new L
      #L.cl <- eVec.cl %*% diag(sqrt(eVal.cl))
      #whichL <- apply(abs(L.cl), 1, which.max) #Which L is each obs in cluster cl assocaiated with?
      #newClust <- whichL
    #}else{
    #  newClust <- rep(1, nrow(A.cl))
    #}

    #return(newClust)
  #})

  #clust2 <- clust
  #for(i in 1:length(allClust)){
  #  clust2[clust == allClust[i]] <- subClusts[[i]]
  #}

  #clust <- interaction(clust, clust2)

  #############
  #SUB-BREAK 2
  #############

  allClust <- unique(clust)

  #Now, for each cluster split based on next SIGN of next largest eigenvector
  subClusts <- lapply(allClust, function(cl){
    #Get cluster specific A
    A.cl <- A[clust == cl, clust == cl, drop = F]

    #As long as cluster isn't singleton or double, try to break down
    if(nrow(A.cl) > 1){

      if(nrow(A.cl) > maxK){
        eDec.cl <- Rfast::eigen.sym(A.cl, k = min(maxK, nrow(A.cl)))
        eVec.cl <- eDec.cl$vectors
        eVal.cl <- eDec.cl$values[,1]
      }else{
        eDec.cl <- eigen(A.cl, symmetric = TRUE)
        eVec.cl <- eDec.cl$vectors
        eVal.cl <- eDec.cl$values
      }
      eVal.cl <- pmax(eVal.cl, rep(0, length(eVal.cl)))
      L.cl <- eVec.cl %*% diag(sqrt(eVal.cl))

    #Get second largest vector
    L2 <- L.cl[,2]
    split <- sign(L2)


      #Evaluate possible split using a version of Frobenius norm if there is an option
      if(!all(diff(split) == 0)){

        #Split if absolute mean of both sides is greater than .5.
        if(mean(abs(L2[split == 1])) > .5 & mean(abs(L2[split == -1])) > .5){
          split <- split
        }else{
          split <- rep(-1, length(L2))
        }

        newClust <- (split + 3)/2
      }else{
        newClust <- rep(1, nrow(A.cl))
      }
    }else{
      newClust <- rep(1, nrow(A.cl))
    }

    return(newClust)
  })

  clust2 <- as.character(clust)
  for(i in 1:length(allClust)){
    clust2[clust == allClust[i]] <- subClusts[[i]]
  }

  clust <- interaction(clust, clust2)

  ###############
  #BUILD - MAIN
  ###############
  clust <- droplevels(clust) #Remove empty levels
  allClust <- names(sort(table(clust), decreasing = TRUE)) #Get unique clusters sorted from largest to smallest

  #Return if all are in the same cluster
  if(length(allClust) == 1){
    return(clust)
  }

  clustPairs <- combn(allClust, 2) #pairwise combinations of all clusters
  clustPairs <- clustPairs[,clustPairs[1,] != clustPairs[2,], drop = F] #No clusters trying to join to themselves

  #Drop pairs which separated in step two above
  pairsWithoutLast <- gsub('.[0-9]$', '', clustPairs)
  clustPairs <- clustPairs[, pairsWithoutLast[1,] != pairsWithoutLast[2,], drop = F]

  #canJoin <- vector('logical', ncol(clustPairs))
  checked <- vector('logical', ncol(clustPairs))

  #While loop for the building
  while(!all(checked)){
    #Identify which pair we're looking at
    p <- which(!checked)[1]
    p1 <- clustPairs[1,p]
    p2 <- clustPairs[2,p]

    #Get points in the two clusters
    cl1 <- clust == p1
    cl2 <- clust == p2

    #Get joining A
    A.cl <- A[cl1|cl2, cl1|cl2, drop = F]

    #Calculate Overlap in 1st two dimensions
    if(nrow(A.cl) > maxK){
      eDec.cl <- Rfast::eigen.sym(A.cl, k = min(maxK, nrow(A.cl)))
      eVec.cl <- eDec.cl$vectors
      eVal.cl <- eDec.cl$values[,1]
    }else{
      eDec.cl <- eigen(A.cl, symmetric = TRUE)
      eVec.cl <- eDec.cl$vectors
      eVal.cl <- eDec.cl$values
    }
    eVal.cl <- pmax(eVal.cl, rep(0, length(eVal.cl)))
    L.cl <- eVec.cl %*% diag(sqrt(eVal.cl))
    L2 <- L.cl[,2]
    split <- sign(L2)

    whichL <- apply(abs(L.cl), 1, which.max)

    #If norm1 < norm2, merge p1 into p2, mark all p1 as checked, mark mustJoin as false
    if(all(whichL == 1) & (mean(L2[split == 1|split == 0]) < .5) & (mean(L2[split == -1 | split == 0]) < .5) ){
      clust[clust == p1] <- p2
      checked[clustPairs[1,] == p1] <- TRUE
    }else{#Otherwise, mark current iteration as checked and add norm1 to joinNorm
      checked[p] <- TRUE
    }

    #Update changing parameters
    clust <- droplevels(clust)
    allClust <- names(sort(table(clust)))
  }

  ##########################
  #BUILD - IF BELOW MINIMUM
  ##########################
  clust <- droplevels(clust) #Remove empty levels
  allClust <- names(sort(table(clust))) #Get unique clusters sorted from smallest to largest

  clustPairs <- combn(allClust, 2) #pairwise combinations of all clusters
  clustPairs <- clustPairs[,clustPairs[1,] != clustPairs[2,], drop = F] #No clusters trying to join to themselves

  #Vectors to save info for build loop
  mustJoin <- sort(table(clust)) < minSize
  #canJoin <- vector('logical', ncol(clustPairs))
  joinNorm <- vector('numeric', ncol(clustPairs))
  checked <- vector('logical', ncol(clustPairs))

  #While loop for the building
  while(any(mustJoin)){
    #Identify which pair we're looking at
    p <- which(!checked)[1]
    p1 <- clustPairs[1,p]
    p2 <- clustPairs[2,p]

    #Get points in the two clusters
    cl1 <- clust == p1
    cl2 <- clust == p2

    #Get joining A
    A.cl <- A[cl1, cl2, drop = F]

    #Calculate Norms
    norm1 <- mean((A.cl - 1)^2)
    norm2 <- mean((A.cl - 0)^2)

    #If norm1 < norm2, merge p1 into p2, mark all p1 as checked, mark mustJoin as false
    if(FALSE){#norm1 < norm2){
      clust[clust == p1] <- p2
      checked[clustPairs[1,] == p1] <- TRUE
      mustJoin[allClust == p1] <- FALSE
    }else{#Otherwise, mark current iteration as checked and add norm1 to joinNorm
      checked[p] <- TRUE
      joinNorm[p] <- norm1
    }

    #If mustJoin still true and all pairs with p1 been checked, merge into best option
    if(mustJoin[allClust == p1] & all(checked[clustPairs[1,] == p1])){
      best <- which.min(joinNorm[clustPairs[1,] == p1])
      clust[clust == p1] <- clustPairs[2,clustPairs[1,] == p1][best]
    }

    #Update changing parameters
    clust <- droplevels(clust)
    allClust <- names(sort(table(clust)))
    mustJoin <- sort(table(clust)) < minSize
  }


  return(clust)
}
