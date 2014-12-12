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
options(error=recover)

library(plyr)
library(phangorn)


#### select 50 samples from full data set

# #read OTU count data
# otu <- t( read.table("../fodor/otu_table_psn_v35.txt", header=T, sep="\t", row.names=1, check.names=FALSE) )

# #read metadata
# id <- read.table("../fodor/v35_map_uniquebyPSN.txt", header=TRUE, sep="\t", row.names=1)

# #examine gut samples only
# body_site <- "Stool"
# site.id <- rownames(id)[which(id$HMPbodysubsite==body_site)]

# #put the samples in columns, and the OTUs in rows
# site <- otu[rownames(otu) %in% site.id,]

# # turn matrix to numeric, transpose so that colnames are sample names, rownames are OTUs
# site.num <- apply(site, 1, function(x){as.numeric(x)})
# rownames(site.num) <- colnames(site)
# site <- site.num

# #pick 50 random samples
# randomSampleIndex <- as.integer(sample(seq(1,length(colnames(site)),1),50,replace=FALSE))
# site.rand <- site[,randomSampleIndex]
# metadata <- id[match(colnames(site.rand),rownames(id)),]

# # remove all 0 count OTUs
# data.sum <- apply(site.rand, 1, sum)
# data <- site.rand[data.sum > 0,]

# #sort from least to most abundant
# data.sum <- apply(data, 1, sum)
# data <- data[order(data.sum),]

# # get rid of extra OTUs in tree
# absent <- tree$tip.label[!(tree$tip.label %in% rownames(data))]
# if (length(absent) != 0) {
# 		tree <- drop.tip(tree, absent)
# }
# write.tree(tree,file="rep_set_v35_subtree.tre")

# #read metadata
# id <- read.table("../fodor/v35_map_uniquebyPSN.txt", header=TRUE, sep="\t", row.names=1)
# id <- id[match(colnames(data),rownames(id)),]
# write.table(id,file="50_random_hmp_gut_samples_metadata.txt",quote=FALSE,sep="\t")

# #save sampled data set
# write.table(data,file="50_random_hmp_gut_samples_otu.txt",quote=FALSE,sep="\t")
# #read with
# # read.table("50_random_hmp_gut_samples_otu.txt", header=TRUE, sep="\t", row.names=1)

####

################# CLR VS PROPORTION TEST ##################
print("clr vs proportion test")

#get proportional data
data.prop <- prop(data)

#get clr data
data.clr <- clr(data)

#create directory for results
mainDir <- getwd()
dir.create(file.path(mainDir, "Compositional_vs_Euclidean"))

summary.prop <- checkMetaData(data.prop,data,metadata,"Compositional_vs_Euclidean/Proportional","unifrac",tree,"proportionsPlots.pdf","Proportions")
summary.clr <- checkMetaData(data.clr,data,metadata,"Compositional_vs_Euclidean/Compositional","euclidean",tree,"CLRPlots.pdf","CLR")


################# SEQUENCING DEPTH ADJUSTMENT TEST ##################
print("sequencing depth adjustment test")

#create directory for results
mainDir <- getwd()
dir.create(file.path(mainDir, "Sequencing_Depth_Adjustments"))

#bootstrap rarification
data.bs <- bootstrap2000(data)
#filter out metadata from discarded samples
metadata.bs <- metadata[match(rownames(data.bs),rownames(metadata)),]
#try proportional and CLR
data.bs.prop <- prop(data.bs)
data.bs.clr <- clr(data.bs)

summary.bs.prop <- checkMetaData(data.bs.prop,data.bs,metadata.bs,"Sequencing_Depth_Adjustments/Bootstrap_Proportions","unifrac",tree,"bootstrapProportionsPlots.pdf","Bootstrap Proportions")

summary.bs.clr <- checkMetaData(data.bs.clr,data.bs,metadata.bs,"Sequencing_Depth_Adjustments/Boostrap_CLR","euclidean",tree,"bootstrapCLRPlots.pdf","Bootstrap CLR")


#jackknife rarification
data.jk <- jackknife2000(data)
#filter out metadata from discarded samples
metadata.jk <- metadata[match(rownames(data.jk),rownames(metadata)),]
#try proportional and CLR
data.jk.prop <- prop(data.jk)
data.jk.clr <- clr(data.jk)

summary.jk.prop <- checkMetaData(data.jk.prop,data.jk,metadata.jk,"Sequencing_Depth_Adjustments/Jackknife_Proportions","unifrac",tree,"jackknifeProportionsPlots.pdf","Jackknife Proportions")

summary.jk.clr <- checkMetaData(data.jk.clr,data.jk,metadata.jk,"Sequencing_Depth_Adjustments/Jackknife_CLR","euclidean",tree,"jackknifeCLRPlots.pdf","Jackknife CLR")


### NEED TO TREAT THIS LIKE REPLICATES -- BUGGY

#dirichlet
data.d <- dirichlet(data)
#try proportional and CLR
data.d.counts <- data.d[[1]]
data.d.prop <- data.d[[2]]
data.d.clr <- clr(data.d.counts)

summary.d.prop <- checkMetaData(data.d.prop,data.d.counts,metadata,"Sequencing_Depth_Adjustments/Dirichlet_Proportions","unifrac",tree,"dirichletProportionsPlots.pdf","Dirichlet Proportions")

summary.d.clr <- checkMetaData(data.d.clr,data.d.counts,metadata,"Sequencing_Depth_Adjustments/Dirichlet_CLR","euclidean",tree,"dirichletCLRPlots.pdf","Dirichlet CLR")


################# SPARSITY TEST ##################
print("sparsity test")

#remove OTUs that are rarer than 1% (1/100) in any every sample
data.sparse1 <- data[,which(apply(data.prop, 2, max) >= 0.01)]
#get proportional data
data.sparse1.prop <- prop(data.sparse1)

#NEED TO TREAT THIS LIKE REPLICATES -- BUGGY

#get clr dirichlet data
data.sparse1.d <- dirichlet(prior(data.sparse1))
data.sparse1.d.prop <- data.sparse1.d[[1]]
data.sparse1.d.counts <- data.sparse1.d[[2]]
data.sparse1.clr <- clr(data.sparse1.d.prop)

summary.sparse1.prop <- checkMetaData(data.sparse1.prop,data.sparse1,metadata,"Sparsity/Less_rare_than_1_percent_proportions","unifrac",tree,"sparse0.1PercentProportionsPlots.pdf","Proportions sparsity filtered at 0.1% abundance")

summary.sparse1.clr <- checkMetaData(data.sparse1.clr,data.sparse1.d.counts,metadata,"Sparsity/Less_rare_than_1_percent_clr","euclidean",tree,"sparse0.1PercentCLRPlots.pdf","CLR sparsity filtered at 0.1% abundance")


#remove OTUs that are rarer than 0.1% (1/1,000) in any every sample
data.sparse01 <- data[apply(data.prop, 1, max) >= 0.001,]
#get proportional data
data.sparse01.prop <- prop(data.sparse01)

#NEED TO TREAT THIS LIKE REPLICATES -- BUGGY

#get clr dirichlet data
#get clr dirichlet data
data.sparse01.d <- dirichlet(prior(data.sparse01))
data.sparse01.d.prop <- data.sparse01.d[[1]]
data.sparse01.d.counts <- data.sparse01.d[[2]]
data.sparse01.clr <- clr(data.sparse01.d.prop)

summary.sparse01.prop <- checkMetaData(data.sparse01.prop,data.sparse01,metadata,"Sparsity/Less_rare_than_01_percent_proportions","unifrac",tree,"sparse0.01PercentProportionsPlots.pdf","Proportions sparsity filtered at 0.01% abundance")

summary.sparse01.clr <- checkMetaData(data.sparse01.clr,data.sparse01.d.counts,metadata,"Sparsity/Less_rare_than_01_percent_clr","euclidean",tree,"sparse0.01PercentCLRPlots.pdf","CLR sparsity filtered at 0.01% abundance")


#remove OTUs that are rarer than 0.001% (1/10,000) in any every sample
data.sparse001 <- data[apply(data.prop, 1, max) >= 0.0001,]
#get proportional data
data.sparse001.prop <- prop(data.sparse001)

#NEED TO TREAT THIS LIKE REPLICATES -- BUGGY

#get clr dirichlet data
#get clr dirichlet data
data.sparse001.d <- dirichlet(prior(data.sparse001))
data.sparse001.d.prop <- data.sparse001.d[[1]]
data.sparse001.d.counts <- data.sparse001.d[[2]]
data.sparse001.clr <- clr(data.sparse001.d.prop)

summary.sparse001.prop <- checkMetaData(data.sparse001.prop,data.sparse001,metadata,"Sparsity/Less_rare_than_001_percent_proportions","unifrac",tree,"sparse0.001PercentProportionsPlots.pdf","Proportions sparsity filtered at 0.001% abundance")

summary.sparse001.clr <- checkMetaData(data.sparse001.clr,data.sparse001.d.counts,metadata,"Sparsity/Less_rare_than_001_percent_clr","euclidean",tree,"sparse0.001PercentCLRPlots.pdf","CLR sparsity filtered at 0.001% abundance")


fakeEffectSummary <- list()

fakeEffectSummary[[1]] <- summary.prop
fakeEffectSummary[[2]] <- summary.clr

fakeEffectSummary[[3]] <- summary.bs.prop
fakeEffectSummary[[4]] <- summary.bs.clr
fakeEffectSummary[[5]] <- summary.jk.prop
fakeEffectSummary[[6]] <- summary.jk.clr
fakeEffectSummary[[7]] <- summary.d.prop
fakeEffectSummary[[8]] <- summary.d.clr

fakeEffectSummary[[9]] <- summary.sparse1.prop
fakeEffectSummary[[10]] <- summary.sparse1.clr
fakeEffectSummary[[11]] <- summary.sparse01.prop
fakeEffectSummary[[12]] <- summary.sparse01.clr
fakeEffectSummary[[13]] <- summary.sparse001.prop
fakeEffectSummary[[14]] <- summary.sparse001.clr


# fileNameList <- list()

# fileNameList[[1]] <- "proportionsPlots.pdf"
# fileNameList[[2]] <- "CLRPlots.pdf"

# fileNameList[[3]] <- "bootstrapProportionsPlots.pdf"
# fileNameList[[4]] <- "bootstrapCLRPlots.pdf"
# fileNameList[[5]] <- "jackknifeProportionsPlots.pdf"
# fileNameList[[6]] <- "jackknifeCLRPlots.pdf"
# fileNameList[[7]] <- "dirichletProportionsPlots.pdf"
# fileNameList[[8]] <- "dirichletCLRPlots.pdf"

# fileNameList[[9]] <- "sparse0.1PercentProportionsPlots.pdf"
# fileNameList[[10]] <- "sparse0.1PercentCLRPlots.pdf"
# fileNameList[[11]] <- "sparse0.01PercentProportionsPlots.pdf"
# fileNameList[[12]] <- "sparse0.01PercentCLRPlots.pdf"
# fileNameList[[13]] <- "sparse0.001PercentProportionsPlots.pdf"
# fileNameList[[14]] <- "sparse0.001PercentCLRPlots.pdf"


# dataTitleList <- list()

# dataTitleList[[1]] <- "Proportions"
# dataTitleList[[2]] <- "CLR"

# dataTitleList[[3]] <- "Bootstrap Proportions"
# dataTitleList[[4]] <- "Bootstrap CLR"
# dataTitleList[[5]] <- "Jackknife Proportions"
# dataTitleList[[6]] <- "Jackknife CLR"
# dataTitleList[[7]] <- "Dirichlet Proportions"
# dataTitleList[[8]] <- "Dirichlet CLR"

# dataTitleList[[9]] <- "Proportions sparsity filtered at 0.1% abundance"
# dataTitleList[[10]] <- "CLR sparsity filtered at 0.1% abundance"
# dataTitleList[[11]] <- "Proportions sparsity filtered at 0.01% abundance"
# dataTitleList[[12]] <- "CLR sparsity filtered at 0.01% abundance"
# dataTitleList[[13]] <- "Proportions sparsity filtered at 0.001% abundance"
# dataTitleList[[14]] <- "CLR sparsity filtered at 0.001% abundance"

# plotAll(fakeEffectSummary,fileNameList,dataTitleList)

save(fakeEffectSummary,file="summary.dat") # load with load(file="summary.dat")