# version 1 (2017/02/05):  first version
# version 2 (2017/02/06):  randomly sample a class if there is a tie

library(DECIPHER)
value <- "Baseline" # name of the run
K <- 8 # k-mer size
S <- function(N) ceiling(32/(1 + exp(-1.445*(log(N) - log(44.33)))))
B <- 1000 # bootstrap replicates
LOOS <- TRUE # perform leave-out-one-sequence?
formula_type <- "OrigSintax"


allkmers <- mkAllStrings(DNA_BASES, K)

# load the training sequences
train <- readDNAStringSet("/home/adithya/RDP/rdpDataset/RDPClassifier_16S_trainsetNo16_rawtrainingdata/trainset16_022016.fa")

# get the full classification (requires "Root;")
classes <- sapply(strsplit(names(train), "Root;", fixed=TRUE),
	`[`,
	2L)
# nothing below here needs to be run for just plotting

# split sequence into k-mers
kmers <- lapply(as.character(train),
	function(x)
		substring(x, 1:(nchar(x) - K + 1), K:nchar(x)))

# index and sort all pure k-mers
kmers <- lapply(kmers, match, table=allkmers)
kmers <- lapply(kmers, function(x) sort(unique(x[!is.na(x)])))

# if performing actual classification, process test sequences
if (LOOS) {
	testkmers <- kmers
} else {
	# load the test sequences
	test <- readDNAStringSet("~/Downloads/SRR5279174.fasta")
	# eventually will need to collapse duplicates here
	# also need to consider sequence orientation (+/-)
	# (check both orientations with a few sequences, pick the best)
	
	# split sequence into k-mers
	testkmers <- lapply(as.character(test),
		function(x)
			substring(x, 1:(nchar(x) - K + 1), K:nchar(x)))
	
	# index and sort all pure k-mers
	testkmers <- lapply(testkmers, match, table=allkmers)
	testkmers <- lapply(testkmers, function(x) sort(unique(x[!is.na(x)])))
}

# ultimately, will need to compute confidence for each class
# but for now only report confidence for the final class label

# count the number of representatives per class
t <- table(classes)
weight <- 1/t[classes] # weigh each class equally
# compute idf weights
counts <- numeric(length(allkmers))
for (i in seq_along(classes)) {
	counts[kmers[[i]]] <- counts[kmers[[i]]] + weight[i]
}
counts <- length(t)/(1 + counts) # Max/(1 + F)
counts <- log(counts) # idf weight is log(Max/(1 + F))
if(formula_type == '+Formula'){ 
	counts[] <- 1
}

# initialize vectors for the output
predicted <- character(length(testkmers))
predicted[] <- NA_character_
confidence <- numeric(length(testkmers))
confidence[] <- NA_real_

pBar <- txtProgressBar(style=3)
for (i in seq_along(testkmers)) {
	if (length(testkmers[[i]])==0)
		next # no k-mers to test
	
	# sample the same k-mers for all genera
	s <- S(length(testkmers[[i]])) # subsample size
	sampling <- matrix(sample(testkmers[[i]], s*B, replace=TRUE), B, s)
	
	# only match unique k-mers to improve speed
	uSampling <- sort(unique(as.vector(sampling)))
	m <- match(sampling, uSampling)
	
	# lookup IDF weights for sampled k-mers
	myweights <- matrix(counts[sampling], B, s)
	
	# record the matches to each genus
	hits <- matrix(0,
		nrow=length(classes),
		ncol=B)
	for (j in seq_along(classes)) {
		matches <- .Call("intMatch",
			uSampling,
			kmers[[j]],
			PACKAGE="DECIPHER")[m]
		matches <- myweights*matches
		# alternatively, the above two commands can be replaced by:
		#matches <- myweights*(sampling %in% kmers[[j]]) # slightly slower
		
		hits[j,] <- rowSums(matches)
	}
	
	if (LOOS)
		hits[i,] <- -1 # eliminate this sequence from matching
	
	# compute confidence from the number of hits per class
	sumOfHits <- numeric(length(classes))
	davg <- mean(rowSums(myweights))
	for (j in seq_len(B)) {
		mymax <- max(hits[,j])
		w <- which(hits[,j]==mymax)
		if (length(w) > 1) {
			selected <- sample(w, 1)
		} else {
			selected <- w
		}
		# our formula:
	
			#sumOfHits[selected] <- sumOfHits[selected] + mymax/davg
		
			# original SINTAX formula:
			sumOfHits[selected] <- sumOfHits[selected] + 1
		
	}	
	
	t <- tapply(sumOfHits, classes, sum)
	w <- which(t==max(t))
	if (length(w) > 1) {
		selected <- sample(w, 1)
	} else {
		selected <- w
	}
	t <- t[selected]
	
	# record the results for this sequence
	predicted[i] <- names(t)
	confidence[i] <- t/B*100 # scale confidence between [0, 100]
	
	setTxtProgressBar(pBar, i/length(testkmers))
	cat(i,' ')
}

# save the results for all sequences
eval(parse(text=paste("p", value, " = predicted", sep="")))
eval(parse(text=paste("c", value, " = confidence", sep="")))
eval(parse(text=paste("save(p", value, ", c", value, ", file='~/Desktop/TEMP/temp", value, ".RData')", sep="")))

####################################
## Plotting
####################################

# create the plot once
plot(NA,
	xlim=c(0, 1),
	ylim=c(0, 1),
	xlab="Fraction of classifiable sequences classified",
	ylab="Classification rate")
# then run the code below once per dataset:

# load the desired results into memory, then:
predicted <- pBaseline
confidence <- cBaseline

t <- table(classes)
weight <- 1/t[classes] # weigh each class equally

singletons <- which(classes %in% names(t[t==1]))
nonsingletons <- which(classes %in% names(t[t > 1]))

thresholds <- seq(0, 100, 0.1)
OC <- MC <- CC <- numeric(length(thresholds))
for (k in seq_along(thresholds)) {
	# OC = (# singletons classified) / (# singletons total)
	OC[k] <- sum(confidence[singletons] >= thresholds[k])/length(singletons)
	
	# MC = (# non-singletons mis-classified) / (# non-singletons classified)
	w <- nonsingletons[which(confidence[nonsingletons] >= thresholds[k])]
	MC[k] <- sum(weight[w[predicted[w] != classes[w]]])/sum(weight[w])
	
	# CC = (# non-singletons classified) / (# non-singletons total)
	CC[k] <- sum(weight[w])/sum(weight[nonsingletons])
}

lines(CC, OC, col="green")
lines(CC, MC, col="green")

legend("topleft",
	legend=c("Original SINTAX", "RDP", "F + IDF log"),
	col=c("gray", "blue", "green"),
	lty=1)
