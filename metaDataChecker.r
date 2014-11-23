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


pairwiseConditionComparator <- function(otu,otucounts,groups,folderName,analysis,tree) {
	
	data <- list()
	
	conditions <- unique(groups)
	data$nConditions <- length(conditions)*(length(conditions)-1)

	if (analysis == "unifrac") {
		data$distMat <- GUniFrac(otucounts,tree,c(1))$unifracs[,,"d_1"]
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
	print("test")
	for (i in 1:(length(conditions)-1)) {
		for(j in (i+1):length(conditions)) {
			print(paste("i",i,"j",j))
			group1indices <- which(groups==conditions[i])
			group2indices <- which(groups==conditions[j])
			data$groups[[index]] <- groups[c(group1indices,group2indices)]
			data$samples[[index]] <- rownames(otu)[c(group1indices,group2indices)]
			data$pcoa[[index]] <- pcoa(data$distMat[which(rownames(data$distMat) %in% data$samples[[index]]),which(colnames(data$distMat) %in% data$samples[[index]])])
			data$wilcoxinRankSum[[index]] <- t(apply(t(otu[c(group1indices,group2indices)]),1,function(x) { as.numeric(wilcox.test(x[1:length(group1indices)], x[length(group1indices):length(c(group1indices,group2indices))])[3]) } ))
			data$wilcoxinRankSumBH[[index]] <- p.adjust(data$wilcoxinRankSum[[index]], method="BH")
			index <- index + 1
		}
	}
	return(data)
}

getSeparation <- function(comparisonSummary,metadata) {
	for (i in 1:length(comparisonSummary)) {
		comparisonSummary[[i]]$nConditions <- levels(factor(metadata[,i]))
		comparisonSummary[[i]]$separation1 <- list()
		comparisonSummary[[i]]$separation2 <- list()
		if (length(comparisonSummary[[i]]$pcoa)>=1) {
			for (j in 1:length(comparisonSummary[[i]]$pcoa)) {
				print(paste("i",i,"j",j,colnames(metadata)[i]))
				comparisonSummary[[i]]$separation1[[j]] <- getPcoaSeparation(comparisonSummary[[i]]$pcoa[[j]],comparisonSummary[[i]]$groups[[j]],1)
				comparisonSummary[[i]]$separation2[[j]] <- getPcoaSeparation(comparisonSummary[[i]]$pcoa[[j]],comparisonSummary[[i]]$groups[[j]],2)
			}
		}
	}
}

#get average PCoA distance on given component
getPcoaSeparation <- function(pcoa,groups,component) {
	conditions <- unique(groups)
	#conditions[1] vs conditions[2] -- expect pairwise condition comparison
	group1indices <- which(groups==conditions[1])
	group2indices <- which(groups==conditions[2])
	distance <- abs(mean(pcoa$vectors[group1indices,component]) - mean(pcoa$vectors[group2indices,component]))
	return(distance)
}

checkMetaData <- function (otu, otucounts, metadata, folderName,analysis,tree)
{

	nSamples <- length(rownames(metadata))

	#add metadata for made up random condition grouping
	imaginaryGrouping <- as.factor(c(rep("1",floor(nSamples/2)),rep("2",(nSamples - floor(nSamples/2)))))
	#shuffle
	imaginaryGrouping <- sample(imaginaryGrouping)
	metadata$imaginaryGrouping <- imaginaryGrouping

	#add metadata for fake "read count category"
	readCountGroups <- as.factor(c(rep("low",floor(nSamples/2)),rep("high",(nSamples - floor(nSamples/2)))))
	readCounts <- apply(data,1,sum)
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
	comparisonData$visitno <- pairwiseConditionComparator(otu,otucounts,metadata$visitno,folderName,analysis,tree)
	comparisonData$sex <- pairwiseConditionComparator(otu,otucounts,metadata$sex,folderName,analysis,tree)
	comparisonData$readCountGroups <- pairwiseConditionComparator(otu,otucounts,metadata$readCountGroups,folderName,analysis,tree)
	comparisonData$imaginaryGrouping <- pairwiseConditionComparator(otu,otucounts,metadata$imaginaryGrouping,folderName,analysis,tree)
	conditionIndices <- c(2,3,9,10)
	comparisonSummary <- getSeparation(comparisonData[conditionIndices],metadata[conditionIndices])

}
