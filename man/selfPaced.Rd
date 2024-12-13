% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/selfPacedCE.R
\name{selfPaced}
\alias{selfPaced}
\title{Self-Paced Clustering Ensemble (SPCE) Algorithm}
\usage{
selfPaced(SList, K, threshold = 0.1)
}
\arguments{
\item{SList}{A list of \verb{n x n} co-association matrices, each representing a base clustering.
Each matrix contains values between 0 and 1 indicating the co-occurrence likelihood of samples.}

\item{K}{An integer specifying the number of clusters to form.}

\item{threshold}{A numeric hyperparameter value (default = 0.1) controlling sparsity internal to the algorithm. We have noticed issues with algorithmic convergence when this value is set too high. If the algorithm is taking super long, try lowering the threshold.}
}
\value{
A vector of cluster assignments for each sample, with length equal to the number of rows (samples) in the input matrices.
Each element is an integer from 1 to \code{K}, indicating the cluster assignment of the corresponding sample.
}
\description{
Implements the SPCE algorithm described in the paper \emph{Self-Paced Clustering Ensemble} by Zhou et al. (2020).
The algorithm integrates multiple base clusterings into a consensus clustering using a self-paced learning strategy.
}
\references{
Zhou, X., Fu, Y., Liang, Z., & Wu, F. (2020). Self-Paced Clustering Ensemble. \emph{IEEE Transactions on Knowledge and Data Engineering}.
}
