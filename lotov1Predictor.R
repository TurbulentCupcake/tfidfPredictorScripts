# In version 1 of loto, we have a balanced training set where
# in our training set. for every prediction, we have to leave out
# the entire taxon. In order to do this, we find out the genus of
# our sequence and check the training set for the indices with that
# genus remove that genus from the entire thing. 


args = (commandArgs(TRUE))

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

load('rdpDataframe.RData')
loadfilename <- paste(c('tfidf',k,'mers.RData'),collapse = "")
load(loadfilename)

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


#  start and end will be split by the number of unique sequences/

# We now proceed to get the unique seqeunces. So in order to do this, we
# get the indices of every unique sequence and take any one of them into our 
# testing sequence (here we take the first one for convenience).

testingSeqsIndices <- lapply(uniqueRank, function(x) { 
		which(x == rank)[1]
	})
testingSeqsIndices <- unlist(testingSeqsIndices)
testingSeqs <- mers[testingSeqsIndices]

# Now that we have our testing sequences, we have to set our training sequences
# to be the set of seqeuences without the genus of our testing seqeunce.

tfidfVals <- eval(parse(text = paste(c('tfidf',k,'mers'), collapse = '')))
predictionVector <- vector(mode = 'character', length = length(testingSeqs))


for(i in start:end) { 

	tfidfSeq <- tfidfVals[[testingSeqsIndices[i]]]
	testSeq <- testingSeqs[i]
	# Now that we have the test genus, we to edit our training sequences such that
	# we remove those sequences that dont have the same sequences as the genus
	testGenus <- uniqueRank[i]
	testGenusIndices <- which(rank %in% testGenus)
	training_db_seqs <- mers[-testGenusIndices]
	training_db_rank <- rank[-testGenusIndices]

	overlapVector <- sapply(training_db_seqs, k = testSeq, FUN = function(X,k) {
			t1 = unlist(k)
			t2 = unlist(X)
	
			# length(intersect(t1,t2))
			sum(tfidfSeq[intersect(t1, t2)])
			})

			maxPos <- which(overlapVector == max(overlapVector))

			if(length(maxPos) > 1) {
					maxPos <- sample(maxPos)[1]
			} else { 
					maxPos <- maxPos[1]
				}


		predicted <- training_db_rank[maxPos]
		predictionVector[i] <- predicted


		cat('Number : ', i ,'\n')
		cat('Test Rank ', testGenus,'\n')
		cat('Predicted ', predicted, '\n')

}

savelink <- paste(c('LOTOv1Prediction_',k,'mers_',end,'_v3.RData'), collapse = "")
save(predictionVector, savelink)