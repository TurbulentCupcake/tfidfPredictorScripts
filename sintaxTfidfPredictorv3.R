

# source(tfidf.R)

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


load('rdpDataframe.RData')
loadfilename <- paste(c('tfidf',k,'mersRdata/tfidf',k,'mers.RData'),collapse = "")

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

	tfidfVals <- eval(parse(text = paste(c('tfidf',k,'mers'), collapse = '')))

	predictionVector <- vector(mode = 'character', length = length(query_ranks))

	for(i in start:end)
	{

		tfidfSeq <- tfidfVals[[i]]
		testSeq <- query_seqs[i]
		testRank <- query_ranks[i]
		training_db_rank <- query_ranks[-i]
		training_db_seqs <- query_seqs[-i]
		testSeq <- unlist(testSeq)
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
				for (i in seq_len(length(w) - 1)) {
					if ((m[w[i + 1]] - m[w[i]]) > (w[i + 1] - w[i])) {
						eliminate[w[i + 1]] <- TRUE
						eliminate[w[i]] <- TRUE
					}
				}

				m[eliminate] = NA_integer_
				matches <- which(!is.na(m) & !eliminate)
				return(length(matches))

			})

			maxPos <- which(overlapVector == max(overlapVector))

			if(length(maxPos) > 1) {
					maxPos <- sample(maxPos)[1]
			} else { 
					maxPos <- maxPos[1]
				}


		predicted <- training_db_rank[maxPos]
		predictionVector[i] <- predicted

		cat('Number : ', i ,'\n')
		cat('Test Rank ', testRank,'\n')
		cat('Predicted ', predicted, '\n')
	}

	savelink <- paste(c('SINTAXtfidf_predictions_',end,'.RData'), collapse = "")
	save(predictionVector, file = savelink)


	
