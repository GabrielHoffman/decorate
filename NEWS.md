# decorate 1.0.31
May 17, 2023
 
 - fix bug `https://github.com/GabrielHoffman/decorate/issues/2`

# decorate 1.0.30
Sept 13, 2022

	- Update dependencies

# decorate 1.0.29
July 23, 2021

	- small updates to documentation with `pkgdown`

# decorate 1.0.28
July 8, 2021

	- incorporate new faster version of adjclust package (>= 0.6.0)

# decorate 1.0.27
November 25, 2019

	- runOrderedClusteringGenome() requires rownames(X) to be found in names(gr)
	- replace cat() with message()

# decorate 1.0.26
September 9, 2019

	- in boxM_permute() return NA's if more variables than samples in the smallest group

# decorate 1.0.25
September 6, 2019

	- use density plot for plotCorrDecay()

# decorate 1.0.24
August 22, 2019

	- add options for plotCorrDecay()

# decorate 1.0.23
August 5, 2019

	- fix bug in evaluateCorrDecay() when chromosome has no peaks

# decorate 1.0.22
August 2, 2019

	- in evalDiffCorr() return effect size estimates
	- in sle.test() and delaneau.test()
		if variable is factor with multiple levels, use Kruskal-Wallis test
		if variable is factor with 2 levels, use Wilcoxon test
		if variable is continuous, cor.test with spearman
	- methods "Cai.max", "Chang.maxBoot", "LC.U", "WL.randProj", "Schott.Frob" now test correlation instead of covariance

# decorate 1.0.21
August 2, 2019

	- fix bug in evalDiffCorr() preventing running "deltaSLE" with multiple categories

# decorate 1.0.20
July 31, 2019

	- modify sle.score() to use lead eigen value of difference matrix
	- add runPermutedData() to get cluster statistics from permuted data
	- faster scoreClusters()
	- plotCorrDecay() 

# decorate 1.0.19
July 26, 2019

	- getClusterRanges() and combineResults() is much faster

# decorate 1.0.18
July 23, 2019

	- add support for the following tests:
	 "Cai.max", "Chang.maxBoot", "LC.U.test", "WL.randProj", "Schott.Frob", "Delaneau", "deltaSLE"
	- add extractCorrelationScores()
	- add getClusterRanges()

# decorate 1.0.17
July 11, 2019

	- reduce memory usage of evaluateCorrDecay()
	- evalDiffCorr() now works on data.frame
	- add corrMatrix.test()
	- add support for following tests: "Box", "Box.permute", "Steiger.fisher", "Steiger", "Jennrich", "Factor", "Mann.Whitney", "Kruskal.Wallis"

# decorate 1.0.16
June 5, 2019

	- improve error message in plotDecorate()

# decorate 1.0.15
May 29, 2019

	- improve error messages for retainClusters()
	- add combineResults()
	- add additional plots to vignette

# decorate 1.0.14
May 28, 2019

	- add plotScatterPairs()
	- plotDecorate() can compute correlation matrix directly from data
	- fixed bug in sign of sLED test statistic
	- plotDecorate() and plotCompareCorr() now have argument absCorr=FALSE as default
	- New simulated data
	- add method "medianCorr" to evalDiffCorr()

# decorate 1.0.13
May 27, 2019

	- apply droplevels() to testVariable in .evalDiffCorr() and plotCompareCorr()

# decorate 1.0.12
March 26, 2019

	- whichCluster() accepts array for feature_id
	- plotDecorate() and plotEnsGenes() now don't print empty page

# decorate 1.0.11
March 25, 2019

	- whichCluster() results now have proper row names
	- add function getFeaturesInClusterList()
	- for plotDecorate(), expand allowable window to 5e7 bp for plotting genes

# decorate 1.0.9
March 22, 2019

	- enable analysys of clusters of size 3

# decorate 1.0.9
March 17, 2019

	- export classes so objects restore correctly from load()

# decorate 1.0.8
March 7, 2019

	- fix error in evaluateCorrDecay()

# decorate 1.0.7
March 4, 2019

	- GRanges object uses names(gr) feature	
	- drop chromosomes from GRanges objects if they don't have any entries

# decorate 1.0.6
Feb 28, 2019

	- fix plotDecorate() for intervals with no genes
	- fix collapseClusters() and plotDecorate() error with feature names
	
# decorate 1.0.5
Feb 27, 2019

	- plotDecorate() now plots genes
	- much faster collapseClusters()

# decorate 1.0.4
Feb 25, 2019

	- Fix small bugs in syntax
	- update vignette to new workflow
	- update readme

# decorate 1.0.3
Feb 22, 2019

	- add retainClusters() to process results of scoreClusters()
	- chance synatax of filterClusters()
	- add collapseCluster() to drop redundant clusters
	- countClusters() fix

# decorate 1.0.2
Feb 21, 2019

	- Add genes to plotDecorate()
	- add handling of multiple clustering cutoff parameters
	- parallelize scoreClusters()

# decorate 1.0.1
	- add evaluateCorrDecay()
	- fix plotting for plotDecorate() and plotCompareCorr()

# decorate 1.0.0
	- Release

# decorate 0.99.3
	- fix issues with ploting trees

# decorate 0.99.2
	- reduce memory usage substantially
	- reduce compute time 100x by using adaptive permutations

# decorate 0.99.1
Jan 22, 2019

    - Initial version
