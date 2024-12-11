#' Robust Spectral Ensemble Clustering
#'
#' @description
#' This function implements Algorithm 1 from 'Robust Spectral Ensemble Clustering via Rank Minimization' by Tao et al. (2019), with some small adjustments found in the code given on Z. Tao's own website at https://ztao.cc/publications.html.
#'
#'
#' @param S Co-association matrix
#' @param K Number of clusters to identify
#' @param lambda1 Tuning parameter 1. From the paper, works best in (.1, 1)
#' @param lambda2 Tuning parameter 2, From the paper, works best in (.1, 1)
#'
#' @return Cluster membership of observations, order matches row/column order in S.
#' @export
robustSpectral <- function(S, K, lambda1, lambda2){
  #Number of observations
  n <- nrow(S)

  #Initialize Convergence Parameters
  rho <- 1.1
  mu <- 1.5/norm(S)
  muM <- 1e10
  tol = 1e-7

  #Initialize Loop parameters
  converged <- FALSE
  q <- 0 #Iteration counter

  #Initialize Matrices
  J <- Z <- E <- Y1 <- Y2 <- Dz <- matrix(0, n, n)
  H <- matrix(0, n, K)
  Lz <- diag(n)
  Dzm12 <- diag(n)

  #Convergence loop
  while(!converged & q < 500){
    #Update J
    J <- svThresholdOp(Z + (1/mu)*Y2,
                       lambda1/mu)

    #Update Z
    part1 <- solve(S %*% t(S) + diag(n))
    part2 <- t(S) %*% S + J - t(S) %*% E
    part3 <- (1/mu)*(t(S) %*% Y1 - Y2 + Dzm12 %*% H %*% t(H) %*% Dzm12)
    Z <- part1 %*% (part2 + part3)

    #Calc Q and column norms
    Q <- S - S %*% Z + Y1/mu
    QcolNorms <- sqrt(colSums(Q^2))

    #Update E
    #THIS IS DIFFERENT FROM THE PAPER, BASED ON ONLINE CODE
    E <- sapply(1:length(QcolNorms), function(qi){
      qn <- QcolNorms[qi]
      Qi <- Q[,qi]
      max((qn - lambda2/mu),0)*Qi/qn
    })

    #Update Dz and Lz
    HZmat <- (Z + t(Z))/2 + H %*% t(H)

    ##################################################################
    #THIS CHUNK IS NOT IN THE PAPER, ADAPTED FROM THE CODE ONLINE ONLY
    #WITHOUT IT THE METHOD FAILS TO CONVERGE
    svd_HZ <- svd(HZmat)
    svd_D <- svd_HZ$d
    r <- sum(svd_D > 1e-4*svd_D[1])
    adj_U <- svd_HZ$u[,1:r] %*% diag(sqrt(svd_D[1:r]))
    adj_U <- t(t(adj_U)/sqrt(rowSums(adj_U^2)))
    HZmat <- (adj_U %*% t(adj_U))^4

    Dz <- diag(rowSums(HZmat))
    Dzm12 <- diag(1/sqrt(diag(Dz)))

    Lz <- Dzm12 %*% HZmat %*% Dzm12

    #Update H
    V_L <- svd(Lz)$u[,1:K]
    H <- Dzm12 %*% V_L
    ##################################################################

    #Calculate convergence criteria (FROM ONLINE CODE)
    crit1 <- (S - S %*% Z) - E
    crit2 <- Z - J
    stopCrit <- max(abs(crit1), abs(crit2))

    converged <- ifelse(stopCrit < tol, T, F)

    #Update Ys
    Y1 <- Y1 + mu*(S - S %*% Z - E)
    Y2 <- Y2 + mu*(Z - J)

    #Update mu
    mu <- min(rho*mu, muM)

    #Update q
    q <- q + 1
  }

  return(kmeans(H, centers = K, iter.max = 30)$cluster)
}
