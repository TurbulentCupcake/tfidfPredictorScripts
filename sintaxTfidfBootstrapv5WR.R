# tfidf contains the tf and idf functions, whose values will be returned 
# into the current namespace. 
# About the functions : 
# ----- tf(query, dataset, mode) 
# This returns the term frequency of the kmers in the query seqeunce as
# a vector, where the names of the vector is the kmers and the values
# of the vector is the term frequency.
# >>>> Possible modes : rawsmooth, lognormal, doublenormK (default k = 0.5) 
# ---- idf(query, dataset, mode) 
# This returns the inverse document frequency as a vector that contains,
# the idf values for every query in our sequence against the dataset.
# The vector names are the names of the kmers and their respective
# idf values. 
# >>>> Possible modes : log, logsmooth, logmax, probfreq
# Refer to tfidf wikipedia page for more information on the tfidf. 

# source('tfidf.R')

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
# loadfilename2 <- paste(c(k,'mersPredictions.RData'), collapse = "")
load(loadfilename)
# load(loadfilename2)
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

	# tfidfVals <- eval(parse(text = paste(c('tfidf',k,'mers'), collapse = '')))
	for(i in start:end)
	{		
		tfidfSeq <- tfidfVals[[i]]
		testSeq <- mers[[i]]
		testRank <- rank[[i]]
		# predictedRankFromPredictions <- predicted_RDP_sintax[i]
		training_db_rank <- rank[-i]
		training_db_seqs <- mers[-i]
		confidenceVector <- vector(mode = 'integer', length = length(uniqueRank))
		names(confidenceVector) <- uniqueRank	
		sequence_df <- data.frame(matrix(NA, nrow = 100,  ncol = 5))
		weights <- tfidfSeq[testSeq]
		probs <- weights/sum(weights)
		samp_matrix_w <- matrix(sample(length(testSeq),3200,replace = TRUE,prob=probs),nrow=100,ncol=32)




		cat('testRank ', testRank, '\n')

		for(j in 1:100) {

			cat('bootstrapNo : ', j, '\n')
			

			testSeq <- unlist(testSeq)
			sampleKmerIndices <- samp_matrix_w[j,]
			# bootstrappedKmers <- testSeq[sampleKmerIndices]
			# The following is our overlap vector, which we can use to find the 
			overlapVector <- sapply(training_db_seqs, k = testSeq, FUN = function(X,k) {
				trainKmers = unlist(X)
				testKmers = unlist(k)

				# replace duplicates with NAs to prevent matching
				trainKmers[trainKmers %in% trainKmers[duplicated(trainKmers)]] <- NA_character_
				testKmers[testKmers %in% testKmers[duplicated(testKmers)]] <- NA_character_

				# find the first and last match
				m <- rep(NA_integer_, length(testKmers))
				m[!is.na(testKmers)] <- match(testKmers[!is.na(testKmers)], trainKmers)


				eliminate <- logical(length(m))
				w <- which(!is.na(m))
				#cat(length(w),'\n')
				if(length(w)!=0) {
					for (i in seq_len(length(w) - 1)) {
						if ((m[w[i + 1]] - m[w[i]]) > (w[i + 1] - w[i])) {
							eliminate[w[i + 1]] <- TRUE
							eliminate[w[i]] <- TRUE
						}
					}
					m[eliminate] = NA_integer_
					matches <- m[sampleKmerIndices]
					matches <- which(!is.na(matches))
					return(length(matches))
				} else { 
					return(0)
				}

						
			})

			# This will return the overlap vector with the hi*di product
			# the next step is to divide them into the 
			maxPos <- which(overlapVector == max(overlapVector))

			cat('Matches with the best = ', max(overlapVector),'\n')

			if(length(maxPos) > 1) {
					maxPos <- sample(maxPos)[1]
			} else { 
					maxPos <- maxPos[1]
				}

			predicted <- training_db_rank[maxPos]
			hi <- max(overlapVector)
			# cat('number of hits = ',hi,'\n')
			cat('Predicted In bootstrap : ', predicted,'\n')
		
		
			# Now that we have the common kmers, we can use the same kmers to get the di from the tfidfSeq	
			# we can use the common kmers to get the values that we need and store it in a dataframe

			# di <- sum(tfidfSeq[bootstrappedKmers])
			sequence_df[j,1] <- predicted
			sequence_df[j,2] <- hi/32
		
			
		}

		
		sequence_df[,5] <- sequence_df[,2]



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
		cat('Final Prediction : ', testRank, '\n')
	}

	# Now we can use the existing bootstrap to begin
	# our predictions using overlapping k-mers

	# Now instead of using only part of the sequence,
	# we will use full length sequences and find out
	# the correct genus using the annotation.

	savelink <- paste(c('confidence_',end,'_v5_WRBalanced.RData'), collapse = "")
	savelink2 <- paste(c('predictions_',end,'_v5_WRBalanced.RData'), collapse = "")
	
	save(bs_confidence_vector, file = savelink)
	save(predictionVector, file = savelink2)

# ----------------------------------------------------------------------------------------
