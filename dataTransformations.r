
#dirichlet method - gets dirichlet distribution of reads per OTU, returns median
# does 128 replicates (median converges at ~64 replicates)
#stolen from Dr. Greg Gloor
rdirichlet <- function (alpha)
{
	test <- data.frame()
	rownames(test) <- c(1,2,3,4,"error")
	n <- 128
  if(length(n) > 1){
  	n <- length(n)
  }
  if(length(n) == 0 || as.integer(n) == 0) {
  	return(numeric(0))
  }
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
	return(t(apply(otu, 1, function(x){x/sum(x)})))
}

#gets clr values from proportional OTU counts, adds 0.5 prior
clr <- function(otu) {
	otu.prior <- prior(otu)
	otu.prop <- prop(otu.prior)
	return(t(apply(otu.prop,1,function(x){log2(x) - mean(log2(x))})))
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
rarefy <- function(otu,minReadCount,withReplacement) {

	#check if OTU names are distinct
	if (length(unique(rownames(otu)))!=length(rownames(otu))) {
		warning("OTU names aren't distinct")
	}


	# put all the OTUs we are sampling from into a list.
	#  OTUs with zero counts aren't included, and OTUs with multiple counts are repeated
	sampleFrom <- list()
	sampleFrom <- lapply(c(1:50),function(x) rep(colnames(otu),otu[x,]))

	#remove reads with less than 2000 counts
	indices <- which(lapply(sampleFrom,function(x) {length(x)})>=2000)
	samples <- sampleFrom[indices]

	#perform bootstrap sampling with replacement
	samples.otu <- lapply(samples,function(x) { sample(x,minReadCount,replace=withReplacement) } )

	#count number of each OTU sampled
	samples.count <- lapply(samples.otu,function(x) {count(x)} )

	#convert OTU names from factors to characters
	bootstrap <- data.frame(lapply(samples.count,function(x) {x$freq[match(colnames(otu),levels(x$x)[as.numeric(x$x)])]} )	)
	bootstrap[is.na(bootstrap)] <- 0

	bootstrap <- t(bootstrap)

	#preserve sample names (so that meta data can be matched) and OTU names
	rownames(bootstrap) <- rownames(otu)[indices]
	colnames(bootstrap) <- colnames(otu)

	#remove all OTUs with zero counts
	bootstrap.sum <- apply(bootstrap,1,sum)
	bootstrap <- bootstrap[which(bootstrap.sum>0),]
	
	return(bootstrap)
}

#gets median of dirichlet distribution for each otu count per sample
dirichlet <- function(otu) {
	return(apply(otu, 2, function(x){rdirichlet(x)}))
}
