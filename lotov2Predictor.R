# In this version of the LOTO, we use the same set of sequences for both training
# and testing, so we dont have to remove entire genera, we only have to finish the 
# 


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


testingSeqsIndices <- lapply(uniqueRank, function(x) { 
        which(x == rank)[1]
    })
testingSeqsIndices <- unlist(testingSeqsIndices)
testingSeqs <- mers[testingSeqsIndices]



# Removing the singleton sequences from our reference database.

tfidfVals <- eval(parse(text = paste(c('tfidf',k,'mers'), collapse = '')))

predictionVector <- vector(mode = 'character', length = length(testingSeqs))

for(i in start:end){ 

        tfidfSeq <- tfidfVals[[testingSeqsIndices[i]]]
        testSeq <- testingSeqs[i]
        testRank <- uniqueRank[i]
        training_db_seqs <- testingSeqs[-i]
        training_db_rank <- uniqueRank[-i]
        overlapVector <- sapply(training_db_seqs, k = testSeq, FUN = function(X,k) {
            t1 = unlist(k)
            t2 = unlist(X)
    
            # length(intersect(t1,t2))
            sum(tfidfSeq[intersect(t1, t2)])
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


savelink <- paste(c('LOTOv2Prediction_',end,'.RData'), collapse = "")
save(predictionVector, file = savelink)


    
