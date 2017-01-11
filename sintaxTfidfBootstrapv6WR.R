

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
		training_db_rank <- rank[-i]
		training_db_seqs <- mers[-i]
		confidenceVector <- vector(mode = 'integer', length = length(uniqueRank))
		names(confidenceVector) <- uniqueRank	
		sequence_df <- data.frame(matrix(NA, nrow = 100,  ncol = 5))
		weights <- tfidfSeq[testSeq]
		probs <- weights/sum(weights)
		samp_matrix_w <- matrix(sample(length(testSeq),3200,replace = TRUE,prob=probs),nrow=100,ncol=32)

		# In the bidirectional sintax version, we do sintax in two ways, first with train and test, then
		# test with train. So essentially. we would have to 




	}





