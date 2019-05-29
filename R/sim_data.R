
# # # http://www.math-evry.cnrs.fr/logiciels/bald
# library(BALD)
# library(GenomicRanges)
# library(mvtnorm)
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


# # Simulate variable to split dataset by
# set.seed(1)
# metadata = data.frame( Disease = sort(factor(sample(0:1, N, replace=TRUE))))


# file = '/Users/gabrielhoffman/workspace/repos/decorate/data/'

# simData_v1 = t(sim$X)
# rownames(simData_v1) = paste('peak', 1:nrow(simData_v1), sep='_')
# colnames(simData_v1) = paste('sample', 1:ncol(simData_v1), sep='_')

# # Simulate first set
# ####################
# # get correlation structure
# C = cor(t(simData_v1))

# # simulate correlation structure with continuous variables
# simData = t(rmvnorm(sum(metadata$Disease == 0), sigma=C))
# rownames(simData) = paste('peak', 1:nrow(simData), sep='_')
# colnames(simData) = paste('sample', 1:ncol(simData), sep='_')


# # Simulate second set
# ####################

# idx = seq(betas$blockSizes[1]+1, betas$blockSizes[1]+betas$blockSizes[2])
# C[idx,idx] = 0
# diag(C) = 1

# simData2 = t(rmvnorm(sum(metadata$Disease == 1), sigma=C))
# rownames(simData2) = paste('peak', 1:nrow(simData2), sep='_')
# colnames(simData2) = paste('sample', 1:ncol(simData2), sep='_')

# # Combine
# #########

# simData = cbind(simData, simData2)

# # C = cor(t(simData))
# # image(C)


# peakIDs = c("peak_75", "peak_76", "peak_77", "peak_78")
# cor(t(simData[peakIDs[1:4], metadata$Disease==0]))
# cor(t(simData[peakIDs[1:4], metadata$Disease==1]))






# # simLocation = GRanges("chr1", IRanges(1:nrow(simData), 1:nrow(simData)), name=rownames(simData))

# # chr20 62044899-62164706
# # custom peak locations
# s = sort(sample(62044899:62164706, nrow(simData), replace=TRUE))
# e = s + c(round(diff(s)/ 2), 10)
# simLocation = GRanges("chr20", IRanges(s,e, names=rownames(simData)))

# # end(simLocation)[-N] > start(simLocation)[-1]

# save(simData, simLocation, metadata, file=paste0(file, "decorateData.rda"))



