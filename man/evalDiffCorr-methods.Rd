% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_sled.R
\docType{methods}
\name{evalDiffCorr}
\alias{evalDiffCorr}
\alias{evalDiffCorr,EList,ANY,GRanges,list-method}
\alias{evalDiffCorr,EList,ANY,GRanges,list,ANY,ANY,ANY,ANY,ANY,ANY,ANY-method}
\alias{evalDiffCorr,matrix,ANY,GRanges,list-method}
\alias{evalDiffCorr,matrix,ANY,GRanges,list,ANY,ANY,ANY,ANY,ANY,ANY,ANY-method}
\alias{evalDiffCorr,data.frame,ANY,GRanges,list-method}
\alias{evalDiffCorr,data.frame,ANY,GRanges,list,ANY,ANY,ANY,ANY,ANY,ANY,ANY-method}
\title{Evaluate Differential Correlation}
\usage{
evalDiffCorr(
  epiSignal,
  testVariable,
  gRanges,
  clustList,
  npermute = c(100, 10000, 1e+05),
  adj.beta = 0,
  rho = 0,
  sumabs.seq = 1,
  BPPARAM = bpparam(),
  method = c("sLED", "Box", "Box.permute", "Steiger.fisher", "Steiger", "Jennrich",
    "Factor", "Mann.Whitney", "Kruskal.Wallis", "Cai.max", "Chang.maxBoot", "LC.U",
    "WL.randProj", "Schott.Frob", "Delaneau", "deltaSLE"),
  method.corr = c("pearson", "kendall", "spearman")
)

\S4method{evalDiffCorr}{EList,ANY,GRanges,list}(
  epiSignal,
  testVariable,
  gRanges,
  clustList,
  npermute = c(100, 10000, 1e+05),
  adj.beta = 0,
  rho = 0,
  sumabs.seq = 1,
  BPPARAM = bpparam(),
  method = c("sLED", "Box", "Box.permute", "Steiger.fisher", "Steiger", "Jennrich",
    "Factor", "Mann.Whitney", "Kruskal.Wallis", "Cai.max", "Chang.maxBoot", "LC.U",
    "WL.randProj", "Schott.Frob", "Delaneau", "deltaSLE"),
  method.corr = c("pearson", "kendall", "spearman")
)

\S4method{evalDiffCorr}{matrix,ANY,GRanges,list}(
  epiSignal,
  testVariable,
  gRanges,
  clustList,
  npermute = c(100, 10000, 1e+05),
  adj.beta = 0,
  rho = 0,
  sumabs.seq = 1,
  BPPARAM = bpparam(),
  method = c("sLED", "Box", "Box.permute", "Steiger.fisher", "Steiger", "Jennrich",
    "Factor", "Mann.Whitney", "Kruskal.Wallis", "Cai.max", "Chang.maxBoot", "LC.U",
    "WL.randProj", "Schott.Frob", "Delaneau", "deltaSLE"),
  method.corr = c("pearson", "kendall", "spearman")
)

\S4method{evalDiffCorr}{data.frame,ANY,GRanges,list}(
  epiSignal,
  testVariable,
  gRanges,
  clustList,
  npermute = c(100, 10000, 1e+05),
  adj.beta = 0,
  rho = 0,
  sumabs.seq = 1,
  BPPARAM = bpparam(),
  method = c("sLED", "Box", "Box.permute", "Steiger.fisher", "Steiger", "Jennrich",
    "Factor", "Mann.Whitney", "Kruskal.Wallis", "Cai.max", "Chang.maxBoot", "LC.U",
    "WL.randProj", "Schott.Frob", "Delaneau", "deltaSLE"),
  method.corr = c("pearson", "kendall", "spearman")
)
}
\arguments{
\item{epiSignal}{matrix or EList of epigentic signal.  Rows are features and columns are samples}

\item{testVariable}{factor indicating two subsets of the samples to compare}

\item{gRanges}{GenomciRanges corresponding to the rows of epiSignal}

\item{clustList}{list of cluster assignments}

\item{npermute}{array of two entries with min and max number of permutations}

\item{adj.beta}{parameter for sLED}

\item{rho}{a large positive constant such that A(X)-A(Y)+diag(rep(rho,p)) is positive definite. Where p is the number of features}

\item{sumabs.seq}{sparsity parameter}

\item{BPPARAM}{parameters for parallel evaluation}

\item{method}{"sLED", "Box", "Box.permute", "Steiger.fisher", "Steiger", "Jennrich", "Factor", "Mann.Whitney", "Kruskal.Wallis", "Cai.max", "Chang.maxBoot", "LC.U", "WL.randProj", "Schott.Frob", "Delaneau", "deltaSLE"}

\item{method.corr}{Specify type of correlation: "pearson", "kendall", "spearman"}
}
\value{
list of result by chromosome and clustList
}
\description{
Evaluate Differential Correlation between two subsets of data
}
\details{
Correlation sturucture between two subsets of the data is evaluated with sparse-Leading-Eigenvalue-Driven (sLED) test:

Zhu, Lingxue, Jing Lei, Bernie Devlin, and Kathryn Roeder. 2017. Testing high-dimensional covariance matrices, with application to detecting schizophrenia risk genes. Annals of Applied Statistics. 11:3 1810-1831. doi:10.1214/17-AOAS1062
}
\examples{
library(GenomicRanges)
library(EnsDb.Hsapiens.v86)

# load data
data('decorateData')

# load gene locations
ensdb = EnsDb.Hsapiens.v86

# Evaluate hierarchical clsutering
treeList = runOrderedClusteringGenome( simData, simLocation ) 

# Choose cutoffs and return clusters
treeListClusters = createClusters( treeList, method = "meanClusterSize", meanClusterSize=c( 10, 20) )

# Plot correlations and clusters in region defined by query
query = range(simLocation)

# Plot clusters
plotDecorate( ensdb, treeList, treeListClusters, simLocation, query)

# Evaluate Differential Correlation between two subsets of data
sledRes = evalDiffCorr( simData, metadata$Disease, simLocation, treeListClusters, npermute=c(20, 200, 2000))

# get summary of results
df = summary( sledRes )

# print results
head(df)

# extract peak ID's from most significant cluster
peakIDs = getFeaturesInCluster( treeListClusters, df$chrom[1], df$cluster[1], "20")

# plot comparison of correlation matrices for peaks in peakIDs
#  where data is subset by metadata$Disease
main = paste0(df$chrom[1], ': cluster ', df$cluster[1])
plotCompareCorr( simData, peakIDs, metadata$Disease) + ggtitle(main)

}
