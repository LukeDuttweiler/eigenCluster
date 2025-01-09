#Least Squares clustering
#'@export
leastSquares <- function(SList){
  #Make co-association matrix
  S <- Reduce('+', SList)/length(SList)

  #Get least squares measure per clustering
  lSqs <- sapply(SList, function(Si){
    return(sum((Si - S)^2))
  })

  #Identify final clustering
  finS <- SList[[which(lSqs == min(lSqs))[1]]]

  g <- igraph::graph_from_adjacency_matrix(finS)
  return(igraph::components(g)$membership)
}
