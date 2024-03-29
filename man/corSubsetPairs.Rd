% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports_parents.R
\name{corSubsetPairs}
\alias{corSubsetPairs}
\title{Compute correlations between pairs of features}
\usage{
corSubsetPairs(
  Y,
  idxi,
  idxj,
  method = c("pearson", "spearman"),
  silent = FALSE,
  setNANtoZero = FALSE
)
}
\arguments{
\item{Y}{matrix where rows are features}

\item{idxi}{indecies}

\item{idxj}{indecies}

\item{method}{specify which correlation method: "pearson" or "spearman"}

\item{silent}{suppress messages}

\item{setNANtoZero}{replace NAN correlation values with a zero}
}
\value{
Compute local correlations between for all k: cor(Y[,idxi[k]], Y[,idxj[k]])
}
\description{
Compute correlations between pairs of features given in idxi and idxj
}
\examples{
# Simulate simple dataset
N = 600
Y = matrix(rnorm(N*100), 100, N)

# select pairs to compute correlations between
i1 = sample.int(N, 200, replace=TRUE)
i2 = sample.int(N, 200, replace=TRUE)

# evaluate all piars
C = corSubsetPairs(t(Y), i1,i2)

# show value
C[i1[10], i2[10]]

# show values from evaluating this pair directly
cor(Y[,i1[10]], Y[,i2[10]])

}
