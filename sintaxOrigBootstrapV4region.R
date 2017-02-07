
load('rdpDataframe.RData')
load('RDP_V4_region.RData')


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
# loadfilename<- paste0(c(k,"mersOrigPredictions.RData"), collapse = "")

rank <- rdp$genus[index]
names(rank) <- rank
sequences <- V4region

uniqueRank <- unique(rank)
names(uniqueRank) <- uniqueRank


mers <- lapply(sequences,
	function (x) {
		substring(x,
			1:(nchar(x) - k + 1L),
			k:nchar(x))
	})
names(mers) <- rank


# --------------------------  THE ALGORTHM -----------------------------------
	

	# Singleton sequences must be present in our query database but 
	# not in our 
 
	

	query_ranks <- rank
	query_seqs <- mers

	# Removing the singleton sequences from our reference database.

	refernece_db_ranks <- rank
	refernece_db_seqs <- mers

	bs_confidence_vector <- vector(mode = 'integer', length=length(rank))
	names(bs_confidence_vector) <- rank
	predictionVector <- vector(mode = 'character', length = length(rank))

	


	for(i in start:end)
	{
		testSeq <- mers[[i]]
		testRank <- rank[[i]]
		# predictedRankFromPredictions <- predictions[i]
		training_db_rank <- rank[-i]
		training_db_seqs <- mers[-i]
		confidenceVector <- vector(mode = 'integer', length = length(uniqueRank))
		names(confidenceVector) <- uniqueRank	
		sequence_df <- data.frame(matrix(NA, nrow = 100,  ncol = 5))


		cat('testRank ', testRank, '\n')


		for(j in 1:100) {

			cat('bootstrapNo : ', j, '\n')
			

			testSeq <- unlist(testSeq)
			sampleKmerIndices <- sample(length(testSeq), s, replace = FALSE)
			bootstrappedKmers <- testSeq[sampleKmerIndices]

			overlapVector <- sapply(training_db_seqs, k = bootstrappedKmers, FUN = function(X,k) {
				t1 = unlist(X)
				t2 = unlist(k)
				length(intersect(t1, t2))
				
			})
			
			overlapVector <- unlist(overlapVector)

			maxPos <- which(overlapVector == max(overlapVector))

			if(length(maxPos) > 1) {
					maxPos <- sample(maxPos)[1]
			} else { 
					maxPos <- maxPos[1]
				}


			predicted <- training_db_rank[maxPos]
			hi <- max(overlapVector)
			cat('Prediction : ', predicted,'\n')


			sequence_df[j,1] <- predicted
			sequence_df[j,2] <- hi/32
			
		}


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

	# Now we can use the existing bootstrap to begin
	# our predictions using overlapping k-mers

	# Now instead of using only part of the sequence,
	# we will use full length sequences and find out
	# the correct genus using the annotation.

	savelink <- paste(c('confidence_modded_',k,'mers_',end,'_V4region.RData'), collapse = "")
	savelink2 <- paste(c('prediction_modded_',k,'mers_',end,'_V4region.RData'), collapse = "")

	save(bs_confidence_vector, file = savelink)
	save(predictionVector, file = savelink2)

# ----------------------------------------------------------------------------------------
