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
	
	data <- list()
	
	conditions <- unique(groups)
	data$nConditions <- length(conditions)*(length(conditions)-1)

	if (analysis == "unifrac") {
		data$distMat <- GUniFrac(t(otu),metaDataChecker.tree,c(1))
	}
	else {
		if (analysis != "euclidean") {
			warning("no valid analysis method specified. valid methods are unifrac and euclidean. proceeding with euclidean")
		}
		data$distMat <- vegdist(t(otu),method="euclidean")
	}

	data$pcoa <- list()

	data$groups <- list()
	data$samples <- list()
	data$wilcoxinRankSum <- list()
	data$wilcoxinRankSumBH <- list()
	index <- 1
	for (i in 1:length(conditions)-1) {
		for(j in i+1:length(conditions)) {
			group1indices <- which(groups==conditions[i])
			group2indices <- which(groups==conditions[j])
			data$groups[[index]] <- groups[c(group1indices,group2indices)
			data$samples[[index]] <- rownames(otu)[c(group1indices,group2indices)]
			data$pcoa <- pcoa(data$distMat[which(rownames(data$distMat) %in% data$samples[[index]]),which(colnames(data$distMat) %in% data$samples[[index]])])
			data$wilcoxinRankSum[[index]] <- t(apply(t(otu[c(group1indices,group2indices)]),1,function(x) { as.numeric(wilcox.test(x[1:length(group1indices)], x[length(group1indices):length(c(group1indices,group2indices))])[3]) } ))
			data$wilcoxinRankSumBH[[index]] <- p.adjust(data$wilcoxinRankSum[[index]], method="BH")
			index <- index + 1
		}
	}

}

getSeparation <- function(comparisonSummary,metadata) {
	for (i in 1:length(comparisonSummary)) {
		comparisonSummary[[i]]$nConditions <- levels(metadata[i])
		comparisonSummary[[i]]$separation1 <- list()
		comparisonSummary[[i]]$separation2 <- list()
		for (j in 1:length(comparisonSummary[[i]]$pcoa)) {
			comparisonSummary[[i]]$separation1[[j]] <- getPcoaSeparation(comparisonSummary[[i]]$pcoa[[j]],comparisonSummary[[i]]$groups[[j]],1)
			comparisonSummary[[i]]$separation2[[j]] <- getPcoaSeparation(comparisonSummary[[i]]$pcoa[[j]],comparisonSummary[[i]]$groups[[j]],2)
		}
	}
}

#get average PCoA distance on given component
getPcoaSeparation <- function(pcoa,groups,component) {
	conditions <- unique(groups)
	#conditions[1] vs conditions[2] -- expect pairwise condition comparison
	group1indices <- which(groups==conditions[1])
	group2indices <- which(groups==conditions[2])
	return average(pcoa[group1indices,component])
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

	comparisonData <- list()

	# compare $visitno, $sex, $HMPbodysubset (multiple sites -- do a pairwise comparison), $readCountGroups, $imaginaryGrouping
	comparisonData$visitno <- pairwiseConditionComparator(otu,metadata$visitno,folderName)
	comparisonData$sex <- pairwiseConditionComparator(otu,metadata$sex,folderName)
	comparisonData$HMPbodysubset <- pairwiseConditionComparator(otu,metadata$HMPbodysubset,folderName)
	comparisonData$readCountGroups <- pairwiseConditionComparator(otu,metadata$readCountGroups,folderName)
	comparisonData$imaginaryGrouping <- pairwiseConditionComparator(otu,metadata$imaginaryGrouping,folderName)

	comparisonSummary <- getSeparation(comparisonData,metadata)

}
