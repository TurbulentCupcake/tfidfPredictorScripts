
load('rdpDataframe.RData')

args = (commandArgs(TRUE))

if(length(args)==0){
    print("No arguments supplied.")
    ##supply default values
    start = 1
    end = 13212
    k = 8    
}else{
    for(i in 1:length(args)){
         eval(parse(text=args[[i]]))
    }
}


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

# Sequences are split into known and unknown sequences, so if your sequence
# has a representative in the rest of the database. So basically there should
# exist another representative of a single genus/ 

# there can exist two types of error, misclassifications and over classifications. 
# If a novel sequence or a sequence that has only one representative is classified into
# some genus, then that is called over classficiation and when a known genus, or
# a genus with other representatives in the reference database, then that is a 
# misclassfiication. 

# These are the indices that will give us our training set

# Controlling for singletons

# All the singletonGenera are supposed tp be novel sequences, so when we take o
# our training set, we should technically have some combination of both singletons and 
# sequneces with more than one. 

# --------------------------  THE ALGORTHM -----------------------------------
	
	# Singleton sequences must be present in our query database but 
	# not in our 
 
	query_ranks <- rank
	query_seqs <- mers


	# Removing the singleton sequences from our reference database.

	# tfidfVals <- tfidf4mers

	predictionVector <- vector(mode = 'character', length = length(query_ranks))

	for(i in start:end)
	{

		# tfidfSeq <- tfidfVals[[i]]
		testSeq <- query_seqs[i]
		testRank <- query_ranks[i]
		training_db_rank <- query_ranks[-i]
		training_db_seqs <- query_seqs[-i]
		testSeq <- unlist(testSeq)
		overlapVector <- sapply(training_db_seqs, k = testSeq, FUN = function(X,k) {
				t1 = unlist(k)
				t2 = unlist(X)
				length(intersect(t1,t2))
		})
		maxPos <- which(overlapVector == max(overlapVector))

		predicted <- training_db_rank[maxPos]
		predictionVector[i] <- predicted

		cat('Number : ', i ,'\n')
		cat('Test Rank ', testRank,'\n')
		cat('Predicted ', predicted, '\n')
	}

	savelink <- paste(c('predictionVector_RDP_sintax_',k,'_',end,'.RData'), collapse = "")
	save(predictionVector, file = savelink)


	
