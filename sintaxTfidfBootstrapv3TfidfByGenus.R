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

# Replcae this with the tfidf for full length sequences.
loadfilename <- paste(c('tfidfByGenusFL.RData'),collapse = "")
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


# MODIFICATION TO EXISTING SINTAX ALGORITHM 

 

# --------------------------  THE ALGORTHM -----------------------------------
	

	# Singleton sequences must be present in our query database but 
	# not in our 


	query_ranks <- rank
	query_seqs <- mers

	# Removing the singleton sequences from our reference database.

	refernece_db_ranks <- rank
	refernece_db_seqs <- mers

	bs_confidence_vector <- vector(mode = 'integer', length=length(mers))
	names(bs_confidence_vector) <- rdp$genus
	tfidfVals <- genusTfidfList

	predictionVector <- vector(mode = 'character', length = length(rank))
	
	for(i in start:end)
	{	
		
		testSeq <- mers[[i]]
		testRank <- rank[[i]]
		
		tfidfSeq <- 
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
			sampleKmerIndices <- samp_matrix_w[j,]#sample(length(testSeq), s, replace = FALSE)
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


		uniquePredictions <- unique(sequence_df[,1])
		cvec <- sapply(uniquePredictions, function(x) { 
				sum(sequence_df[which(sequence_df[,1] == x),5])
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

	savelink <- paste(c('confidence_',end,'_v3_RDP.RData'), collapse = "")
	savelink2 <- paste(c('predictions_',end,'_v3_RDP.RData'), collapse = "")

	save(bs_confidence_vector, file = savelink)
	save(predictionVector, file = savelink2)

# ----------------------------------------------------------------------------------------
