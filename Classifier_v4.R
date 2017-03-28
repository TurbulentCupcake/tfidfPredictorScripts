# version 1 (2017/02/05):  first version
# version 2 (2017/02/06):  randomly sample a class if there is a tie
# version 3 (2017/02/07):  select the best representative per class
# version 4 (2017/03/08):  added LOOT option, more parameters, and USEROOT

args = (commandArgs(TRUE))
if(length(args)==0){
    print("No arguments supplied.")
    ##supply default values
    VALUE <- "100v1" # name of the run
    B <- 100 # bootstrap replicates
    LOOS <- TRUE # perform leave-out-one-sequence?
    LOOT <- FALSE # perform leave-out-one-taxon?
    SINTAX <- FALSE # no IDF weights, S=f(N), or formula
    TOPHIT <- FALSE # only use the top hit from each genus
	USEROOTS <- FALSE # re-balance by rooting 

}else{
    for(i in 1:length(args)){
         eval(parse(text=args[[i]]))
    }
}

library(DECIPHER)
options(warn=3) # error if a warning occurs
# will eventually want to set.seed() for repetability

TRAINFILE <- "trainset16_022016.fa"
# VALUE <- "100v1" # name of the run
VALUE <- paste(c(B,'v',version), collapse = '')
K <- 8 # k-mer size
B <- 100 # bootstrap replicates
LOOS <- TRUE # perform leave-out-one-sequence?
LOOT <- FALSE # perform leave-out-one-taxon?
if (!LOOS && !LOOT)
	TESTFILE <- "~/Downloads/SRR5279174.fasta"
SINTAX <- FALSE # no IDF weights, S=f(N), or formula
if (SINTAX) {	
	S <- function(N) return(32)
} else {
	S <- function(N) ceiling(32/(1 + exp(-1.445*(log(N) - log(44.33)))))
}
# imbalance parameters
TOPHIT <- FALSE # only use the top hit from each genus
USEROOTS <- FALSE # re-balance by rooting 
if (SINTAX && TOPHIT)
	warning("SINTAX and TOPHIT used in conjunction.")
if (LOOS && LOOT)
	stop("Only LOOS or LOOT may be TRUE.")
# Balance sequence by picking one sequence for each genus as the triaing
# set
BALANCE <- FALSE 

allkmers <- mkAllStrings(DNA_BASES, K)

# load the training sequences
train <- readDNAStringSet(TRAINFILE)

if(BALANCE) { 
	seqnames <-  sapply(strsplit(names(train), "Root;", fixed=TRUE),
	`[`,
	2L)
	genera <- unique(seqnames)
	indices <- sapply(genera,FUN = function(x) { 
			which(seqnames == x)[1]
		})
	if(!all(seqnames[indices] == genera)) stop("Error occured here")
	train <- train[indices]

}

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
if (LOOS || LOOT) {
	testkmers <- kmers
} else {
	# load the test sequences
	test <- readDNAStringSet(TESTFILE)
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
classTable <- table(classes)
weight <- 1/classTable[classes] # weigh each class equally
# compute idf weights
counts <- numeric(length(allkmers))
for (i in seq_along(classes)) {
	counts[kmers[[i]]] <- counts[kmers[[i]]] + weight[i]
}
counts <- length(classTable)/(1 + counts) # Max/(1 + F)
counts <- log(counts) # idf weight is log(Max/(1 + F))
if (SINTAX)
	counts[] <- 1

# find roots for each class
if (USEROOTS) {
	f <- function(N) (1/(1 + exp(-1*(log(N) - log(1)))))
	roots <- f(classTable)
} else {
	roots <- rep(1, length(classTable))
}

# obtain a list with the index for each group
nonsingletons <- which(classTable > 1)
groups <- lapply(names(classTable),
	function(x) {
		which(classes==x)
	})

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
	
	if (LOOS) {
		hits[i,] <- -1 # eliminate this sequence from matching
	} else if (LOOT) {
		w <- which(classes==classes[i])
		hits[w,] <- -1 # eliminate this taxon from matching
	}
	
	if (TOPHIT) {
		# find the top hit in each nonsingleton group
		sumHits <- rowSums(hits)
		topHits <- groups # initialize with singletons
		for (j in seq_along(nonsingletons)) {
			w <- which.max(sumHits[groups[[nonsingletons[j]]]])
			topHits[[nonsingletons[j]]] <- topHits[[nonsingletons[j]]][w]
		}
		topHits <- unlist(topHits)
	} else {
		topHits <- seq_along(classes) # use every sequence
	}
	
	# compute confidence from the number of hits per class
	totHits <- numeric(length(topHits))
	davg <- mean(rowSums(myweights))
	for (j in seq_len(B)) {
		mymax <- max(hits[topHits, j])
		w <- which(hits[topHits, j]==mymax)
		if (length(w) > 1) {
			selected <- sample(w, 1)
		} else {
			selected <- w
		}
		
		if (SINTAX) { # use original SINTAX formula
			totHits[selected] <- totHits[selected] + 1
		} else { # use our formula
			totHits[selected] <- totHits[selected] + mymax/davg
		}
	}
	
	if (!TOPHIT) {
		t <- tapply(totHits, classes, sum)
		totHits <- t[names(classTable)] # order the same as classTable
	}
	
	w <- which(totHits==max(totHits))
	if (length(w) > 1) {
		selected <- sample(w, 1)
	} else {
		selected <- w
	}
	
	# record the results for this sequence
	predicted[i] <- names(classTable[selected])
	confidence[i] <- ((totHits[selected]/B)^roots[selected])*100 # scale confidence between [0, 100]
	
	setTxtProgressBar(pBar, i/length(testkmers))
}

# save the results for all sequences
eval(parse(text=paste("p", VALUE, " = predicted", sep="")))
eval(parse(text=paste("c", VALUE, " = confidence", sep="")))
eval(parse(text=paste("save(p", VALUE, ", c", VALUE, ", file='", VALUE, ".RData')", sep="")))

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
predicted <- pBaseline16S
confidence <- cBaseline16S

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

lines(CC, OC, col="black")
lines(CC, MC, col="black")

legend("topleft",
	legend=c("Original SINTAX", "+ Formula (F)", "F + IDF log"),
	col=c("gray", "black", "green"),
	lty=1)
