spectralClust <- function(S, K){
  eigS <- Rfast::eigen.sym(S, k = K)$vectors

  return(kmeans(eigS, centers = K)$cluster)
}
