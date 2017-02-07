

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

loadfilename <- paste(c('VRtfidf',k,'mers.RData'),collapse = "")
load(loadfilename)
load('rdpDataframe.RData')
load('RDP_V4_region.RData')

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


bs_confidence_vector <- vector(mode = 'integer', length=length(rank))
names(bs_confidence_vector) <- uniqueRank
tfidfVals <- eval(parse(text = paste(c('vr',k,'mers'), collapse = '')))




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
		training_db_rank <- rank[-i]
		training_db_seqs <- mers[-i]
		confidenceVector <- vector(mode = 'integer', length = length(uniqueRank))
		names(confidenceVector) <- uniqueRank	
		sequence_df <- data.frame(matrix(NA, nrow = 100,  ncol = 5))
		weights <- tfidfSeq[testSeq]
		probs <- weights/sum(weights)
		samp_matrix_w <- matrix(sample(length(testSeq),3200,replace = TRUE,prob=probs),nrow=100,ncol=32)

		# In the bidirectional sintax version, we do sintax in two ways, first with train and test, then
		# test with train. 
		# So in our first step, we proceed as we did with regular sintax, so we have to bootstrap our test
		# sequence 50 times against our training set as we normally do and assign the confidences
		# as we normally would by dividing the bootstrap replicate by 32.

		# When we do it the other way round, we take our training set as our set of test sequences and we train
		# against our previous test sequence. So what you do is : 
		# for every bootstrap 
		# go through each of sequence in the training set
		# bootstrap 32 out of it 
		# using those 32 try to find out how many hits you get against the test sequeence
		# switch to the next one, take a new bootstrap replicate, match
		# and continue to do this for every single sequence in the 13211 and find the 
		# genus that we hit best on 

		cat('testRank = ', testRank,'\n')

		# FRONT DIRECTION ---------------------------------------------------------------------

		for(j in 1:50) { 

			cat('bootstrapNo : ', j, '\n')


			testSeq <- unlist(testSeq)
			sampleKmerIndices <- samp_matrix_w[j,]
			bootstrappedKmers <- testSeq[sampleKmerIndices]
			# The following is our overlap vector, which we can use to find the 
			overlapVector <- sapply(training_db_seqs, k = bootstrappedKmers, FUN = function(X,k) {
				t1 = unlist(X)
				t2 = unlist(k)
				# So instead of just getting the overlaps as 1s and 0s
				# you would have to find the overlap in the tfidfseq
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
			cat('Predicted In bootstrap : ', predicted,'\n')
		
		
			# Now that we have the common kmers, we can use the same kmers to get the di from the tfidfSeq	
			# we can use the common kmers to get the values that we need and store it in a dataframe

			di <- sum(tfidfSeq[bootstrappedKmers])
			sequence_df[j,1] <- predicted
			sequence_df[j,2] <- hi/32
		

		}

		cat('SWITHCING DIRECTIONS','\n')
		# REVERSE DIRECTION -------------------------------------------------------------

		# We have two versions in this case, either we can take the random sampling of 32 
		# or we can bias according to the kmer. 

		# Bias with kmer ---------------------------------------------------------------

		# we can bias the 32 picked kmers but that process is really slow to complete in
		# reaasonable amount of time.
		for(j in 51:100){

			cat('bootstrapNo : ', j, '\n')


			overlapVector <- sapply(training_db_seqs, k = testSeq, function(X,k) { 
					t1 <- unlist(X)
					t2 <- unlist(k)
					length(intersect(t2,t1[sample(length(t1),32,replace = TRUE)])) 
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
			
			# confidenceVector[predicted] <- confidenceVector[predicted] + 1
			cat('Predicted In bootstrap : ', predicted,'\n')
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


	savelink <- paste(c('confidence_',end,'_v6_VarReg.RData'), collapse = "")
	savelink2 <- paste(c('predictions_',end,'_v6_VarReg.RData'), collapse = "")
	
	save(bs_confidence_vector, file = savelink)
	save(predictionVector, file = savelink2)





