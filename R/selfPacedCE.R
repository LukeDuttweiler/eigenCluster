#' Self-Paced Clustering Ensemble (SPCE) Algorithm
#'
#' Implements the SPCE algorithm described in the paper *Self-Paced Clustering Ensemble* by Zhou et al. (2020).
#' The algorithm integrates multiple base clusterings into a consensus clustering using a self-paced learning strategy.
#'
#' @param SList A list of `n x n` co-association matrices, each representing a base clustering.
#'              Each matrix contains values between 0 and 1 indicating the co-occurrence likelihood of samples.
#' @param K An integer specifying the number of clusters to form.
#' @param threshold A numeric hyperparameter value (default = 0.5) controlling the sparsity of the output. Larger values enforce stricter sparsity.
#'
#' @return A vector of cluster assignments for each sample, with length equal to the number of rows (samples) in the input matrices.
#'         Each element is an integer from 1 to `K`, indicating the cluster assignment of the corresponding sample.
#'
#'
#' @references Zhou, X., Fu, Y., Liang, Z., & Wu, F. (2020). Self-Paced Clustering Ensemble. *IEEE Transactions on Knowledge and Data Engineering*.
#'
#' @export
selfPaced <- function(SList, K, threshold = .5){
  #number of clusterings
  m <- length(SList)
  #number of samples
  n <- nrow(SList[[1]])

  #First main calc
  S.hat <- Reduce('+', SList)/m
  omega <- ((S.hat == 1) | (S.hat == 0)) + 0

  #Initialize Parameters
  S <- S.hat
  D <- colSums(S)
  L <- diag(D) - S
  eigL <- eigen(L, symmetric = TRUE)
  Y <- eigL$vectors[,(n - K + 1):n, drop = F]
  rho <- 1
  alpha <- rep(1/m, m)
  gamma <- threshold^2*m^2

  #Loop over r (speed of learning)
  for(r in c(.9,.8,.7,.6,.5)){
    #Initialize Lambda
    lambda <- 2*(((r-1)^2)*r + (r^2)*(1-r))*m^2

    #Initialize B
    ASq <- lapply(SList, function(Si){
      (S-Si)^2
    })
    B <- Reduce('+', mapply(`/`, ASq, alpha, SIMPLIFY = FALSE))

    #Initialize W
    W <- pmin(lambda/(2*B), 1)

    #iteration indicator
    q <- 1
    while(q < 500){
      #Break if rank is already correct
      rankL <- qr(L)$rank
      if(rankL == n-K){
        break
      }

      #Calculate C
      pq <- expand.grid(1:n, 1:n)
      yNorms <- as.matrix(dist(Y, diag = TRUE))^2
      Sweight <- Reduce('+', mapply(`/`, SList, alpha, SIMPLIFY = FALSE))
      C <- (Sweight - (rho*yNorms/(2*W)))/sum(1/alpha)

      #Calculate tau
      tau <- gamma/Reduce('+', lapply(alpha, function(alph){W^2/alph}))

      #Update S
      S <- C
      S[C < sqrt(tau)] <- 0
      S[C >=1] <- 1

      #Update Y
      D <- colSums(S)
      L <- diag(D) - S
      eigL <- eigen(L, symmetric = TRUE)
      Y <- eigL$vectors[,(n - K + 1):n, drop = F]

      #Update alpha
      d <- sapply(SList, function(Si){
        norm((S - Si)*W, type = 'F')^2
      })
      alpha <- sqrt(d)/sum(sqrt(d))

      #Get rank of L
      rankL <- qr(L)$rank
      if(rankL > n-K){
        rho <- 2*rho
        q <- q+1
      }else if(rankL < n-K){
        rho <- rho/2
        q <- q+1
      }
    }
  }
  S <- (S != 0) + 0

  g <- igraph::graph_from_adjacency_matrix(S)
  return(igraph::components(g))
}
