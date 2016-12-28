# THis is the accompanying code to bootstrap predictions
# for the predictor. 
# 


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


loadfilename <- paste(c('tfidf',k,'mers.RData'),collapse = "")
loadfilename2 <- paste(c('LOTOv1_',k,'mersPredictions.RData'), collapse = "")
load(loadfilename)
load(loadfilename2)
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


testingSeqsIndices <- lapply(uniqueRank, function(x) { 
		which(x == rank)[1]
	})
testingSeqsIndices <- unlist(testingSeqsIndices)
testingSeqs <- mers[testingSeqsIndices]



bs_confidence_vector <- vector(mode = 'integer', length=length(testingSeqs))
names(bs_confidence_vector) <- uniqueRank
tfidfVals <- eval(parse(text = paste(c('tfidf',k,'mers'), collapse = '')))


for(i in start:end) { 

	tfidfSeq <- tfidfVals[[testingSeqsIndices[i]]]
	testSeq <- testingSeqs[i]
	# Now that we have the test genus, we to edit our training sequences such that
	# we remove those sequences that dont have the same sequences as the genus
	testGenus <- uniqueRank[i]
	testGenusIndices <- which(rank %in% testGenus)
	training_db_seqs <- mers[-testGenusIndices]
	training_db_rank <- rank[-testGenusIndices]
	sequence_df <- data.frame(matrix(NA, nrow = 100,  ncol = 5))


	cat('testRank ', testRank, '\n')

		for(j in 1:100) {

			cat('bootstrapNo : ', j, '\n')
			

			testSeq <- unlist(testSeq)
			sampleKmerIndices <- sample(length(testSeq), s, replace = FALSE)
			bootstrappedKmers <- testSeq[sampleKmerIndices]
			# The following is our overlap vector, which we can use to find the 
			overlapVector <- sapply(training_db_seqs, k = bootstrappedKmers, FUN = function(X,k) {
				t1 = unlist(X)
				t2 = unlist(k)
				# So instead of just getting the overlaps as 1s and 0s
				# you would have to find the overlap in the tfidfseq
				sum(tfidfSeq[intersect(t1, t2)])
 					
			})

			# This will return the overlap vector with the hi*di product
			# the next step is to divide them into the 
			maxPos <- which(overlapVector == max(overlapVector))

			if(length(maxPos) > 1) {
					maxPos <- sample(maxPos)[1]
			} else { 
					maxPos <- maxPos[1]
				}

			predicted <- training_db_rank[maxPos]
			hi <- max(overlapVector)
			
			# confidenceVector[predicted] <- confidenceVector[predicted] + 1
			cat('Predicted In bootstrap : ', predicted,'\n')
		
		
			# Now that we have the common kmers, we can use the same kmers to get the di from the tfidfSeq	
			# we can use the common kmers to get the values that we need and store it in a dataframe

			di <- sum(tfidfSeq[bootstrappedKmers])
			sequence_df[j,1] <- predicted
			sequence_df[j,2] <- hi
			sequence_df[j,3] <- di
			
		}


		# Once we have the values for all the bootstraps. we need to calculate te di/davg fore every bootstrap
		# we do this by doing the following 
		davg <- sum(sequence_df[,3])/100
		sequence_df[,4] <- sequence_df[,3]/davg

		# In this version, our confidence is determined by the fact that we tie the
		# goodness of our bootstrap with how good our prediction was in that bootstrap

		sequence_df[,5] <- sequence_df[,4] * (sequence_df[,2]/sequence_df[,3])


		# we run into an issue where our actual rank may not be present in one of the predicted ranks,
		# if thats the case, then we need to control for it by assigning the value of tha bootstrap for that
		# to 0.

		predictedRanks <- sequence_df[,1]
		predictedRanks <- unique(predictedRanks)

		if(predictions[i] %in% predictedRanks) { 
				# getting the threshold for our current sequence.
				confidence <- sum(sequence_df[which(sequence_df[,1] == predictions[i]),5])
			} else if((predictions[i] %in% predictedRanks) == FALSE) { 
					confidence <- 0 
					}

		bs_confidence_vector[i] <- confidence
		cat('Query Seq : ', i,'\n')
		cat('Final Prediction : ', testRank, '\n')
	}


savelink <- paste(c('LOTOv1Bootstrap_',k,'mers_',end,'_v3.RData'), collapse = "")
	



