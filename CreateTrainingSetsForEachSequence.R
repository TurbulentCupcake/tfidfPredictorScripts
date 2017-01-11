

args = (commandArgs(TRUE))
# Bi-directional sintax
# use 50 replicates for bootstrapping on each size
# use intersect. 



if(length(args)==0){
    print("No arguments supplied.")
    ##supply default values
    start = 1
    end = 13212
    k = 8 	
    s = 32
}else{
    for(i in 1:length(args)){
         eval(parse(text=args[[i]]))
    }
}

loadfilename <- paste(c('tfidf',k,'mers.RData'),collapse = "")
load(loadfilename)
load('rdpDataframe.RData')

rank <- rdp$genus
names(rank) <- rank
sequences <- rdp$sequences

uniqueRank <- unique(rank)
names(uniqueRank) <- uniqueRank


mers <- lapply(sequences,
	function (x) {
		substring(x,
			1:(nchar(x) - k + 1L),
			k:nchar(x))
	})
names(mers) <- rank

# In this program, we have to create a training set for each test sequence in our dataset.
# This means, for every test sequence in our dataset, we will have to create a new smaller 
# dataset containing 2472 for 2472 genera in the dataset where each index corresponds to a 
# sequence from a unique genus from the dataset. 

# How do we do this?
# First, we take a test seqeuence from the dataset and make the rest of the dataset as our
# training set. Then, we segregate the training set based on the genus they represent.
# We can do this by picking out the corresponding indices for each genus.
# Once we have segregated, we train genus-by-genus. So for every genus that is represented
# by our training set, we train against the seqeunces in that genus and find the sequence
# that we hit most on that genus's set of sequence and pick that index and
# add it to our vector of indices for that sequence. 

# Use the v4 sintax

tfidfVals <- eval(parse(text = paste(c('tfidf',k,'mers'), collapse = '')))

findBestIndices <- function(indices, test, train, samp_matrix_w) {

    # we now have the training sequence from each genus
    # Our objective at this point is to do 100 boostrap replicates
    # on the new training set and find the index with the most number of
    # hits. 


    # make a vector keep counts for which is the best vector
    counts <- integer(length(train))

    test <- unlist(test)

    for(j in 1:100) {

     #   cat('Bootstrap No', j,'\n')

         bootstrappedKmers <- test[samp_matrix_w[j,]]
         overlapVector <- sapply(train, k = bootstrappedKmers, FUN = function(x,k) { 
                t1 <- unlist(x)
                t2 <- unlist(k)
                length(intersect(t1,t2))
            })

         maxPos <- which(overlapVector == max(overlapVector))
         if(length(maxPos) > 1){
            maxPos <- sample(maxPos)[1]
         }
         counts[maxPos] <- counts[maxPos] + 1
       #  cat('counts', counts)
    }
    names(counts) <- indices

    maxIndices <- names(counts[which(counts == max(counts))])
    if(length(maxIndices) > 1){
        maxIndices <- sample(maxIndices)[1]
    }

    cat('predictedIndex = ',maxIndices,'\n')
    return(maxIndices)
}

indexHolder <- as.list(integer(13212))


    for(i in start:end) { 


        # Pick the test sequence

        tfidfWeights <- tfidfVals[[i]]
        testSeq <- mers[[i]]
        testRank <- rank[[i]]
        trainingSeqs <- mers[-i]
        trainingRank <- rank[-i]
        indexWithNames <- seq(1,13211)
        names(indexWithNames) <- trainingRank
        weights <- tfidfWeights[testSeq]
        probs <- weights/sum(weights)
        samp_matrix_w <- matrix(sample(length(testSeq),3200,replace = TRUE,prob=probs),nrow=100,ncol=32)

        cat('Sequence No : ', i,'\n')
        cat('Sequence Genus = ',testRank,'\n')


        # Get the indices of seqeunces of every genus

        uniqueTrainingRanks <- unique(trainingRank)
        segIndex <- lapply(uniqueTrainingRanks, function(x) { 
                which(x == names(indexWithNames))
            })

        # Find the best index from processing every genus

        genusIndices <- sapply(segIndex, FUN = function(x) { 

                indices <- unlist(x)
                train <- trainingSeqs[indices]
                test <- unlist(testSeq)
                returnIndex <-  findBestIndices(indices, test, train, samp_matrix_w)
                cat('Returned Index = ', returnIndex,'\n')
                return(returnIndex)
            })

        genusIndices <- strtoi(genusIndices)
        genusIndices[which(genusIndices > i)] <- genusIndices[which(genusIndices > i)] + 1
        names(genusIndices) <- uniqueTrainingRanks
        indexHolder[[i]] <- genusIndices

        # cat('Final Indices', genusIndices)

    }