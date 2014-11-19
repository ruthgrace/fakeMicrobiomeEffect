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
site.rand <- site[,as.integer(sample(seq(1,length(colnames(site)),1),50,replace=FALSE))]

# remove all 0 count OTUs
data.sum <- apply(site.rand, 1, sum)
data <- site.rand[data.sum > 0,]

#save sampled data set
write.table(data,file="50_random_hmp_gut_samples_otu.txt",quote=FALSE,sep="\t")
#read with
# read.table("50_random_hmp_gut_samples_otu.txt", header=TRUE, sep="\t", row.names=1)

################# CLR VS PROPORTION TEST ##################

#get clr data

#get proportional data
data.prop <- apply(data, 2, function(x){x/sum(x)})

################# SPARSITY TEST ##################



	# #remove OTUs that are rarer than 0.1% in any every sample
	# site.abund <- site.non0[apply(site.prop0, 1, max) >= 0.001,]
	# site.prop <- apply(site.abund, 2, function(x){x/sum(x)})

