#' Lancichinetti and Fortunato Consensus Clustering using base algorithms from iGraph
#'
#' This function performs L-F consensus clustering on a weighted adjacency matrix using
#' a specified community detection algorithm (default is 'louvain'). The function iterates
#' until a consensus partition is achieved or a maximum number of iterations is reached.
#'
#' Based on ``Consensus clustering in complex networks'', Lancichinetti & Fortunato, 2012.
#'
#' @param adjMat A square, symmetric adjacency matrix representing the network, with edge weights between nodes.
#' @param np Integer. The number of partitions to generate in each iteration (default is 20).
#' @param threshold Numeric. A threshold for filtering edges; values below this are set to zero (default is 0.2).
#' @param alg Character. The community detection algorithm to use, with 'louvain' as the default.
#'             The algorithm should be compatible with the `igraph` package, ie. `cluster_`alg.
#' @param maxIter Integer. The maximum number of iterations to achieve consensus (default is 100).
#' @param ... Additional arguments passed to the community detection algorithm specified in `alg`.
#'
#' @return A vector indicating the final partition for each node, with nodes in the same community
#'         assigned the same value.
#' @details The function iteratively applies the clustering algorithm and adjusts the adjacency matrix
#'          until consensus is reached across partitions or `maxIter` iterations are completed.
#' @export
lf.consensus <- function(adjMat,
                         np = 20,
                         threshold = .2,
                         alg = 'louvain',
                         maxIter = 100,
                         ...){
  #adjMat[adjMat < threshold] <- 0

  #set up indicators
  t <- 1
  complete <- FALSE

  #Loop until all partitions agree or until max iter is reached
  while(t <= maxIter & !complete){
    adjMat <- lf.oneIter(adjMat = adjMat,
                         np = np,
                         threshold = threshold,
                         alg = alg,
                         ... = ...)

    #If all iterations agree, move on
    if(all(adjMat == 1| adjMat == 0)){
      complete <- TRUE
    }else{
      adjMat[adjMat < threshold] <- 0
      t <- t+1
    }
  }

  #Make into partition vector from matrix
  o <- do.call('order', as.data.frame(adjMat))
  adjMatO <- adjMat[o,]
  partO <- cumsum(!duplicated(adjMatO))
  partition <- partO[order(o)]
  print(paste0('Number of Iterations: ', t))
  return(partition)
}

#' Single Iteration of iGraph clustering to support LF Consensus Clustering
#'
#' Helper function to perform a single iteration of consensus clustering on an adjacency matrix
#' by generating multiple partitions and averaging them into a new adjacency matrix.
#'
#' @inheritParams lf.consensus
#'
#' @return A new adjacency matrix, updated based on averaged partitions.
#' @details The function creates multiple partitions, adjusts singleton nodes to improve connectivity,
#'          and computes an averaged adjacency matrix across partitions.
lf.oneIter <- function(adjMat,
                       np,
                       threshold,
                       alg,
                       ...){
  adjGraph <- igraph::graph_from_adjacency_matrix(adjMat, mode = 'undirected', weighted = TRUE, diag = F)

  clustAlg <- eval(parse(text = paste0('igraph::cluster_', alg)))

  partitions <- lapply(1:np, function(i){
    origClust <- clustAlg(adjGraph, ...)$membership

    singles <- which(table(origClust) == 1)

    for(s in singles){
      k <- which(origClust == s)
      #Skip if we've added things together already
      if(length(k) > 1){
        next
      }
      #If no skip, add to the closest option
      origClust[k] <- origClust[which.max(adjMat[k,-k])]
    }
    return(as.factor(origClust))
  })

  newMats <- lapply(partitions, function(p){
    return(outer(p,p, FUN = '==')+0)
  })

  newMat <- Reduce('+', newMats)/np

  return(newMat)
}
