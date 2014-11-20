##############################################################
# 
# Ruth Grace Wong (ruthgracewong@gmail.com)
# Gloor Lab, University of Western Ontario
#
# function checks metadata to see what yields the best separation
# for a fake effect
#
# outputs plots & pcoa values into folder,
# with summary of variance explained
#
###############################################################

library(GUniFrac)
library(ape)
library(phangorn)

metaDataChecker.init <- function(treeFilePath) {
	metaDataChecker.tree <<- read.tree(treeFilePath)
	if (!is.rooted(tree)) {
		metaDataChecker.tree <<- midpoint(metaDataChecker.tree)
	}
}

pairwiseConditionComparator <- function(otu,groups,folderName,analysis) {
	
	otu <- list()
	
	conditions <- unique(groups)
	otu$nConditions <- length(conditions)*(length(conditions)-1)

	if (analysis == "unifrac") {
		otu$distMat <- GUniFrac(t(otu),metaDataChecker.tree,c(1))
	}
	else {
		if (analysis != "euclidean") {
			warning("no valid analysis method specified. valid methods are unifrac and euclidean. proceeding with euclidean")
		}
		otu$distMat <- vegdist(t(otu),method="euclidean")
	}

	otu$pcoa <- pcoa(otu$distMat)

	otu$groups <- list()
	otu$samples <- list()
	index <- 1
	for (i in 1:length(conditions)-1) {
		for(j in i+1:length(conditions)) {
			otu$groups[[index]] <- 
			otu$samples[[index]] <- 
			index <- index + 1
		}
	}	

	# add wilcoxin rank sum

}

getSeparation <- function(comparisonSummary,groups) {
	for (i in 1:length(comparisonSummary)) {
		comparisonSummary[[i]]$nConditions <- length(unique(groups))
		for (j in 1:length(comparisonSummary[[i]]$pcoa)) {
			comparisonSummary[[i]]$separation1 <- getPcoaSeparation(comparisonSummary[[i]]$pcoa[[j]],comparisonSummary[[i]]$groups[[j]],1)
			comparisonSummary[[i]]$separation2 <- getPcoaSeparation(comparisonSummary[[i]]$pcoa[[j]],comparisonSummary[[i]]$groups[[j]],2)
		}
	}
}

#get average PCoA distance on first vector
getPcoaSeparation <- function(pcoa,groups,component) {
	conditions <- unique(groups)

	# not done yet
}

checkMetaData <- function (otu, metadata, folderName)
{

	nSamples <- length(rownames(metadata))

	#add metadata for made up random condition grouping
	imaginaryGrouping <- as.factor(c(rep("1",floor(nSamples/2)),rep("2",(nSamples - floor(nSamples/2)))))
	#shuffle
	imaginaryGrouping <- sample(imaginaryGrouping)
	metadata$imaginaryGrouping <- imaginaryGrouping

	#add metadata for fake "read count category"
	readCountGroups <- as.factor(c(rep("low",floor(nSamples/2)),rep("high",(nSamples - floor(nSamples/2)))))
	readCounts <- apply(data,2,sum)
	readCountOrder <- order(readCounts)
	#assign the half of samples with lowest counts to "low" read count condition
	readCountGroups[readCountOrder[c(1:floor(nSamples/2))]] <- levels(readCountGroups)[2]
	#assign the half of samples with highest counts to "high" read count condition
	readCountGroups[readCountOrder[c((floor(nSamples/2)+1):nSamples)]] <- levels(readCountGroups)[1]
	metadata$readCountGroups <- readCountGroups

	#make output folder if it doesn't exist already
	mainDir <- getwd()
	dir.create(file.path(mainDir, folderName))

	comparisonSummary <- list()

	# compare $visitno, $sex, $HMPbodysubset (multiple sites -- do a pairwise comparison), $readCountGroups, $imaginaryGrouping
	comparisonSummary$visitno <- pairwiseConditionComparator(otu,metadata$visitno,folderName)
	comparisonSummary$sex <- pairwiseConditionComparator(otu,metadata$sex,folderName)
	comparisonSummary$HMPbodysubset <- pairwiseConditionComparator(otu,metadata$HMPbodysubset,folderName)
	comparisonSummary$readCountGroups <- pairwiseConditionComparator(otu,metadata$readCountGroups,folderName)
	comparisonSummary$imaginaryGrouping <- pairwiseConditionComparator(otu,metadata$imaginaryGrouping,folderName)

	comparisonSummary <- getSeparation(comparisonSummary)

}
