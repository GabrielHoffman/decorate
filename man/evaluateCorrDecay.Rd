% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evaluateCorrDecay.R
\name{evaluateCorrDecay}
\alias{evaluateCorrDecay}
\title{Evaluate the decay of correlation versus distance between features}
\usage{
evaluateCorrDecay(treeList, gr, chromArray = seqlevels(gr), verbose = TRUE)
}
\arguments{
\item{treeList}{list of hclust objects}

\item{gr}{GenomicRanges object corresponding to features clustered in treeList}

\item{chromArray}{Use this only this set of chromosmes.  Can substantially reduce memory usage}

\item{verbose}{show progress}
}
\value{
a data.frame of distance and correlation value for all pairs of features already evalauted in treeList.  Note that runOrderedClusteringGenome() that returns treeList only evalutes correlation between a specified number of adjacent peaks
}
\description{
For pairs of features evaluate the physical distance and the correlation
}
\examples{
library(GenomicRanges)
library(ggplot2)

data('decorateData')

# Evaluate hierarchical clustering
treeList = runOrderedClusteringGenome( simData, simLocation ) 

# Evaluate how correlation between features decays with distance
dfDist = evaluateCorrDecay( treeList, simLocation )

# make plot
plotCorrDecay( dfDist )

}
