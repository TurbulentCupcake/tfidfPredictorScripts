
load('rdpDataframe.RData')

args = (commandArgs(TRUE))

if(length(args)==0){
    print("No arguments supplied.")
    ##supply default values
    start = 1
    end = 20000
    k = 8
    s = 32
}else{
    for(i in 1:length(args)){
         eval(parse(text=args[[i]]))
    }
}


loadfilename<-paste0(c(k,"mersPredictionsOrig.RData"), collapse = "")
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

	


	for(i in start:end)
	{
		testSeq <- mers[i]
		testRank <- rank[i]
		predictedRankFromPredictions <- predictions[i]
		training_db_rank <- rank[-i]
		training_db_seqs <- mers[-i]
		confidenceVector <- vector(mode = 'integer', length = length(uniqueRank))
		names(confidenceVector) <- uniqueRank	

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
			
			predicted <- training_db_rank[which(overlapVector == max(overlapVector))[1]]
			confidenceVector[predicted] <- confidenceVector[predicted] + 1
			cat('Predicted In bootstrap : ', predicted,'\n')
		}

		bs_confidence_vector[i] <- confidenceVector[predictedRankFromPredictions]
		cat('Query Seq : ', i,'\n')
		cat('Final Prediction : ', testRank, '\n')
	}

	# Now we can use the existing bootstrap to begin
	# our predictions using overlapping k-mers

	# Now instead of using only part of the sequence,
	# we will use full length sequences and find out
	# the correct genus using the annotation.

	savelink <- paste(c('origSintaxBootstrap_',k,'mers_s',s,'_',end,'.RData'), collapse = "")
	
	save(bs_confidence_vector, file = savelink)

# ----------------------------------------------------------------------------------------
