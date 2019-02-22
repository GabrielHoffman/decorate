
# # http://www.math-evry.cnrs.fr/logiciels/bald
# library(BALD)
# library(GenomicRanges)
# set.seed(5)
# N <- 1000
# r2 <- 0.4
# corr <- 0.5
# blockSizes <- sample(5:80, 10, replace=TRUE)
# p <- sum(blockSizes)
# p
# sig.blocks <- 20
# l <- length(sig.blocks)
# if (sum(is.element(sig.blocks,blockSizes))!=l) stop("sig.blocks must be a subset of blockSizes !")
# nb.per.block <- 1
# betas <- simBeta(blockSizes, sig.blocks, nb.per.block)
# ## effective blockSizes used for simulation !
# betas$blockSizes
# betaSNP <- betas$betaMat[, "betaSNP"]
# sim <- simulation(N, betaSNP, betas$blockSizes, corr, r2)


# file = '/Users/gabrielhoffman/workspace/repos/decorate/data/'

# simData = t(sim$X)
# rownames(simData) = paste('peak', 1:nrow(simData), sep='_')
# colnames(simData) = paste('sample', 1:ncol(simData), sep='_')

# # simLocation = GRanges("chr1", IRanges(1:nrow(simData), 1:nrow(simData)), name=rownames(simData))

# # chr20 62044899-62164706
# # custom peak locations
# s = sort(sample(62044899:62164706, nrow(simData), replace=TRUE))
# e = s + c(round(diff(s)/ 2), 10)
# simLocation = GRanges("chr20", IRanges(s,e), name=rownames(simData))

# # end(simLocation)[-N] > start(simLocation)[-1]

# # Simulate variable to split dataset by
# set.seed(1)
# metadata = data.frame( Disease = factor(sample(0:1, ncol(simData), replace=TRUE)))

# save(simData, simLocation, metadata, file=paste0(file, "decorateData.rda"))



