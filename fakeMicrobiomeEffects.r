##############################################################
# 
# Ruth Grace Wong (ruthgracewong@gmail.com)
# Gloor Lab, University of Western Ontario
#
# tries to create fake effect using metadata and different
#  types of analysis of 50 randomly selected healthy gut
#  samples from the Human Microbiome Project.
#
# data is from http://www.hmpdacc.org/HMQCP/
#
###############################################################

library(plyr)

############ HELPER FUNCTIONS ############

#dirichlet method - gets dirichlet distribution of reads per OTU, returns median
#stolen from Dr. Greg Gloor
rdirichlet <- function (n, alpha)
{
  if(length(n) > 1) n <- length(n)
  if(length(n) == 0 || as.integer(n) == 0) return(numeric(0))
  n <- as.integer(n)
  if(n < 0) stop("integer(n) can not be negative in rtriang")

  if(is.vector(alpha)) alpha <- t(alpha)
  l <- dim(alpha)[2]
  x <- matrix(rgamma(l * n, t(alpha)), ncol = l, byrow=TRUE)  # Gere le recycling
  return(x / rowSums(x))
}

#add prior of 0.5 to raw otu counts
prior <- function(otu) {
	otu.prior <- otu
	otu.prior[otu.prior==0] <- 0.5
	return(otu.prior)
}

#gets proportions from OTU counts
prop <- function(otu) {
	return(apply(otu, 2, function(x){x/sum(x)}))
}

#gets clr values from proportional OTU counts, adds 0.5 prior
clr <- function(otu) {
	otu.prior <- prior(otu)
	otu.prop <- prop(otu.prior)
	return(apply(otu.prop,2,function(x){log2(x) - mean(log2(x))}))
}


bootstrap2000 <- function(otu) {
	return(rarefy(otu,2000,TRUE))
}

jackknife2000 <- function(otu) {
	return(rarefy(otu,2000,FALSE))
}

#input OTU table must have OTU names as row names, these must be distinct
#  samples with less than minReadCount reads are discarded,
#  and everything else is bootstrapped/jackknifed to 2000
#  (depending on the withReplacement boolean value)
rarefy <- function(otu,minReadCount,withReplacement)

	#check if OTU names are distinct
	if (length(unique(rownames(otu)))!=length(rownames(otu))) {
		warning("OTU names aren't distinct")
	}

	sampleNames <- colnames(otu)

	# put all the OTUs we are sampling from into a list.
	#  OTUs with zero counts aren't included, and OTUs with multiple counts are repeated
	sampleFrom <- list()
	sampleFrom <- lapply(c(1:50),function(x) rep(rownames(data),data[,x]))

	#remove reads with less than 2000 counts
	indices <- which(lapply(sampleFrom,function(x) {length(x)})>=2000)
	samples <- sampleFrom[indices]

	#perform bootstrap sampling with replacement
	samples.otu <- lapply(samples,function(x) { sample(x,minReadCount,replace=withReplacement) } )

	#count number of each OTU sampled
	samples.count <- lapply(samples.otu,function(x) {count(x)} )

	#convert OTU names from factors to characters
	bootstrap <- data.frame(lapply(samples.count,function(x) {x$freq[match(rownames(otu),as.character(x$x))]} )	)
	bootstrap[is.na(bootstrap)] <- 0

	#preserve sample names (so that meta data can be matched) and OTU names
	rownames(bootstrap) <- rownames(otu)
	colnames(bootstrap) <- sampleNames[indices]

	#remove all OTUs with zero counts
	bootstrap.sum <- apply(bootstrap,1,sum)
	bootstrap <- bootstrap[which(bootstrap.sum>0),]

	return(bootstrap)
}

#gets median of dirichlet distribution for each otu count per sample
dirichlet <- function(otu) {
	return(apply(otu, 2, function(x){rdirichlet(1,x)}))
}

####################### INTITIALIZATION SCRIPT #########################

#read metadata
id <- read.table("../fodor/v35_map_uniquebyPSN.txt", header=TRUE, sep="\t", row.names=1)

#read OTU count data
otu <- t( read.table("../fodor/otu_table_psn_v35.txt", header=T, sep="\t", row.names=1, check.names=FALSE) )

#examine gut samples only
body_site <- "Stool"
site.id <- rownames(id)[which(id$HMPbodysubsite==body_site)]

#put the samples in columns, and the OTUs in rows
site <- otu[rownames(otu) %in% site.id,]

# turn matrix to numeric, transpose so that colnames are sample names, rownames are OTUs
site.num <- apply(site, 1, function(x){as.numeric(x)})
rownames(site.num) <- colnames(site)
site <- site.num

#pick 50 random samples
randomSampleIndex <- as.integer(sample(seq(1,length(colnames(site)),1),50,replace=FALSE))
site.rand <- site[,randomSampleIndex]
metadata <- id[match(colnames(data),rownames(id)),]

# remove all 0 count OTUs
data.sum <- apply(site.rand, 1, sum)
data <- site.rand[data.sum > 0,]

#sort from least to most abundant
data.sum <- apply(data, 1, sum)
data <- data[order(data.sum)]

#save sampled data set
write.table(data,file="50_random_hmp_gut_samples_otu.txt",quote=FALSE,sep="\t")
#read with
# read.table("50_random_hmp_gut_samples_otu.txt", header=TRUE, sep="\t", row.names=1)

fakeEffectSummary <- list()
summaryIndex <- 1

################# CLR VS PROPORTION TEST ##################

#get proportional data
data.prop <- prop(data)

#get clr data
data.clr <- clr(data)

#create directory for results
mainDir <- getwd()
dir.create(file.path(mainDir, "Compositional_vs_Euclidean"))

fakeEffectSummary[[summaryIndex]] <- checkMetaData(data.prop,metadata,"Compositional_vs_Euclidean/Raw_Counts")
summaryIndex <- summaryIndex+1

fakeEffectSummary[[summaryIndex]] <- checkMetaData(data.prop,metadata,"Compositional_vs_Euclidean/Proportional")
summaryIndex <- summaryIndex+1

fakeEffectSummary[[summaryIndex]] <- checkMetaData(data.prop.metadata,"Compositional_vs_Euclidean/Compositional")
summaryIndex <- summaryIndex+1

################# SEQUENCING DEPTH ADJUSTMENT TEST ##################

#create directory for results
mainDir <- getwd()
dir.create(file.path(mainDir, "Sequencing_Depth_Adjustments"))

#bootstrap rarification
data.bs <- bootstrap2000(data)
#try proportional and CLR
data.bs.prop <- prop(data.bs)
data.bs.clr <- clr(data.bs)

[summaryIndex]] <- checkMetaData(data.bs.prop,metadata,"Sequencing_Depth_Adjustments/Bootstrap_Proportions")
summaryIndex <- summaryIndex+1

[summaryIndex]] <- checkMetaData(data.bs.clr,metadata,"Sequencing_Depth_Adjustments/Boostrap_CLR")
summaryIndex <- summaryIndex+1

#jackknife rarification
data.jk <- jackknife(data)
#try proportional and CLR
data.jk.prop <- prop(data.jk)
data.jk.clr <- clr(data.jk)

[summaryIndex]] <- checkMetaData(data.jk.prop,metadata,"Sequencing_Depth_Adjustments/Jackknife_Proportions")
summaryIndex <- summaryIndex+1

[summaryIndex]] <- checkMetaData(data.jk.clr,metadata,"Sequencing_Depth_Adjustments/Jackknife_CLR")
summaryIndex <- summaryIndex+1

#dirichlet
data.d <- dirichlet(data)
#try proportional and CLR
data.d.prop <- prop(data.d)
data.d.clr <- clr(data.d)

[summaryIndex]] <- checkMetaData(data.d.prop,metadata,"Sequencing_Depth_Adjustments/Dirichlet_Proportions")
summaryIndex <- summaryIndex+1

[summaryIndex]] <- checkMetaData(data.d.clr,metadata,"Sequencing_Depth_Adjustments/Dirichlet_CLR")
summaryIndex <- summaryIndex+1


################# SPARSITY TEST ##################

#remove OTUs that are rarer than 1% (1/100) in any every sample
data.sparse1 <- data[apply(data.prop0, 1, max) >= 0.01,]
#get proportional data
data.sparse1.prop <- prop(data.sparse1)
#get clr dirichlet data
data.sparse1.d <- dirichlet(data.sparse1)
data.sparse1.clr <- clr(data.sparse1.d)

#remove OTUs that are rarer than 0.1% (1/1,000) in any every sample
data.sparse01 <- data[apply(data.prop0, 1, max) >= 0.001,]
#get proportional data
data.sparse01.prop <- prop(data.sparse01)
#get clr dirichlet data
data.sparse01.d <- dirichlet(data.sparse01)
data.sparse01.clr <- clr(data.sparse01.d)

#remove OTUs that are rarer than 0.001% (1/10,000) in any every sample
data.sparse001 <- data[apply(data.prop0, 1, max) >= 0.0001,]
#get proportional data
data.sparse001.prop <- prop(data.sparse001)
#get clr dirichlet data
data.sparse001.d <- dirichlet(data.sparse001)
data.sparse001.clr <- clr(data.sparse001.d)




