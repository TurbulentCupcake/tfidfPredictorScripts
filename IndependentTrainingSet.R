

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
load(loadfilename)
load('rdpDataframe.RData')
load('trainingSetForEachSequence.RData')


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


bs_confidence_vector <- vector(mode = 'integer', length=length(rank))
names(bs_confidence_vector) <- uniqueRank
tfidfVals <- eval(parse(text = paste(c('tfidf',k,'mers'), collapse = '')))


# MODIFICATION TO EXISTING SINTAX ALGORITHM 
 

# --------------------------  THE ALGORTHM -----------------------------------
	

	# Singleton sequences must be present in our query database but 
	# not in our 


	query_ranks <- rank
	query_seqs <- mers


	# Removing the singleton sequences from our reference database.

	refernece_db_ranks <- rank
	refernece_db_seqs <- mers

	# bs_confidence_vector <- vector(mode = 'integer', length=length(mers))
	predictionVector <- vector(mode = 'character', length = length(rank))
	#  names(bs_confidence_vector) <- rdp$genus

	for(i in start:end) { 



		tfidfSeq <- tfidfVals[[i]]
		testSeq <- mers[[i]]
		testRank <- rank[[i]]
		# predictedRankFromPredictions <- predicted_RDP_sintax[i]
		training_db_rank <- rank[tsList[[i]]]
		training_db_seqs <- mers[tsList[[i]]]
		confidenceVector <- vector(mode = 'integer', length = length(uniqueRank))
		names(confidenceVector) <- uniqueRank	
		sequence_df <- data.frame(matrix(NA, nrow = 100,  ncol = 5))
		weights <- tfidfSeq[testSeq]
		probs <- weights/sum(weights)
		samp_matrix_w <- matrix(sample(length(testSeq),3200,replace = TRUE,prob=probs),nrow=100,ncol=32)


		for(j in 1:100) { 

					cat('bootstrapNo : ', j, '\n')


			testSeq <- unlist(testSeq)
			sampleKmerIndices <- samp_matrix_w[j,]
			bootstrappedKmers <- testSeq[sampleKmerIndices]
			overlapVector <- sapply(training_db_seqs, k = bootstrappedKmers, FUN = function(X,k) {
				t1 = unlist(X)
				t2 = unlist(k)
				length(intersect(t1, t2))
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
			cat('Match Kmer Count = ', hi,'\n')
			cat('Predicted In bootstrap : ', predicted,'\n')
		
		
			# Now that we have the common kmers, we can use the same kmers to get the di from the tfidfSeq	
			# we can use the common kmers to get the values that we need and store it in a dataframe

			# di <- sum(tfidfSeq[bootstrappedKmers])
			sequence_df[j,1] <- predicted
			sequence_df[j,2] <- hi/32
	
		}


		
		sequence_df[,5] <- sequence_df[,2]
		uniquePredictions <- unique(sequence_df[,1])
		cvec <- sapply(uniquePredictions, function(x) { 
				sum(sequence_df[which(sequence_df[,1] == x),2])
			})
		names(cvec) <- uniquePredictions
		maxPos2 <- which(cvec == max(cvec))
		if(length(maxPos2) > 1) { 
			maxPos2 <- sample(maxPos2)[1]
		}else{
			maxPos2 <- maxPos2[1]
		}

		confidence <- cvec[maxPos2]
		prediction <- uniquePredictions[maxPos2]


		bs_confidence_vector[i] <- confidence
		predictionVector[i] <- prediction
		cat('Query Seq : ', i,'\n')
		cat('Final Prediction : ', testRank, '\n')

	}



	savelink <- paste(c('confidence_',end,'_IndependentTraining.RData'), collapse = "")
	savelink2 <- paste(c('predictions_',end,'_IndependentTraining.RData'), collapse = "")
	
	save(bs_confidence_vector, file = savelink)
	save(predictionVector, file = savelink2)



