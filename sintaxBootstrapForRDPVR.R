
load('rdpDataframe.RData')
load('RDP_V4_region.RData')


args = (commandArgs(TRUE))

if(length(args)==0){
    print("No arguments supplied.")
    ##supply default values
    start = 1
    end = 13212
    k=8
    s=32
}else{
    for(i in 1:length(args)){
         eval(parse(text=args[[i]]))
    }
}

rank <- rdp$genus[index]
names(rank) <- rank
sequences <- rdp$sequences[index]


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

	bs_confidence_vector <- vector(mode = 'integer', length=length(mers))
	names(bs_confidence_vector) <- rdp$genus[index]
	predictionVector <- vector(mode = 'character', length = length(mers))




	for(i in start:end)
	{
		testSeq <- mers[[i]]
		testRank <- rank[[i]]
		# predictedRankFromPredictions <- predicted_RDP_sintax[i]
		training_db_rank <- rank[-i]
		training_db_seqs <- mers[-i]
		confidenceVector <- vector(mode = 'integer', length = length(uniqueRank))
		names(confidenceVector) <- uniqueRank
		sequence_df <- data.frame(matrix(NA, nrow = 100,  ncol = 5))
		samp_matrix_w <- matrix(sample(length(testSeq),3200,replace = TRUE),nrow=100,ncol=32)


		cat('testRank ', testRank, '\n')


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

			# predicted <- training_db_rank[which(overlapVector == max(overlapVector))[1]]
			# confidenceVector[predicted] <- confidenceVector[predicted] + 1
			# cat('Predicted In bootstrap : ', predicted,'\n')

			maxPos <- which(overlapVector == max(overlapVector))

			if(length(maxPos) > 1) {
					maxPos <- sample(maxPos)[1]
			} else {
					maxPos <- maxPos[1]
				}

			predicted <- training_db_rank[maxPos]
			hi <- max(overlapVector)
			cat('Value in bootstrap : ', hi,'\n')

			# confidenceVector[predicted] <- confidenceVector[predicted] + 1
			cat('Predicted In bootstrap : ', predicted,'\n')


			# Now that we have the common kmers, we can use the same kmers to get the di from the tfidfSeq
			# we can use the common kmers to get the values that we need and store it in a dataframe

			sequence_df[j,1] <- predicted
			sequence_df[j,2] <- 1

		}

		# bs_confidence_vector[i] <- confidenceVector[predictedRankFromPredictions]
		# cat('Query Seq : ', i,'\n')
		# cat('Final Prediction : ', testRank, '\n')
		sequence_df[,5] <- sequence_df[,2]


		# we run into an issue where our actual rank may not be present in one of the predicted ranks,
		# if thats the case, then we need to control for it by assigning the value of tha bootstrap for that
		# to 0.

		# prediction <- sample(names(which(table(sequence_df[,1]) == max(table(sequence_df[,1])))))[1]
		# confidence <- table(sequence_df[,1])[prediction]
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
		cat('Actual Sequence : ', testRank,'\n')
		cat('Final Prediction : ', prediction, '\n')
		cat('Final Confidence : ', confidence, '\n')

	}

	# Now we can use the existing bootstrap to begin
	# our predictions using overlapping k-mers

	# Now instead of using only part of the sequence,
	# we will use full length sequences and find out
	# the correct genus using the annotation.

	savelink <- paste(c('sintax_',end,'_confidence_VR.RData'), collapse = "")
	savelink2 <- paste(c('sintax_',end,'_prediction_VR.RData'), collapse = "")


	save(bs_confidence_vector, file = savelink)
	save(predictionVector, file = savelink2)

# ----------------------------------------------------------------------------------------
