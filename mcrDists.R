load('distanceMatrix.RData')

s <- 32
version <- 3
veks <- c('8mers')

for(lab2 in 1:1) { 

loadfile1 <- paste(c("../PredictionsRData/",veks[lab2],"Confidencesv",version,"s",s,".RData"), collapse = "")
loadfile2 <- paste(c("../PredictionsRData/",veks[lab2],"Predictions.RData"), collapse = "")
load(loadfile1)
load(loadfile2)
linetype = lab2
# load('origConfidence.RData')
# load('origPredictions.RData')
load('rdpDataframe.RData')



thresholdValues <- c(0.80)



predictions <- predictions
actual <- rdp$genus
singletonGenera <- names(which(table(actual) == 1))
knownGenera <- names(which(table(actual) != 1))
lengthOfKnown <- length(actual)-length(singletonGenera)
MCvector <- vector(mode = 'double', length = length(thresholdValues))
OCvector <- vector(mode = 'double', length = length(thresholdValues))
bootstrapValues <- confidences
# bootstrapValues <- bootstrapValues/100
unclassifiedVector <- vector(mode = 'double', length = length(thresholdValues))

OCpredictions <- vector(mode = "character", length = 13212)
MCpredictions <- vector(mode = "character", length = 13212)


for(i in seq_along(thresholdValues)) {

	 cat('Threshold = ', thresholdValues[i],'\n')
	 currentThreshold <- thresholdValues[i]

	 #Creating our new vector here 
	 newPredictions <- vector(mode = 'character', length = length(predictions))
	 newPredictions[] <- NA_character_

	 for(j in seq_along(bootstrapValues)) {
	 	if(is.na(bootstrapValues[j]) == FALSE) {
	 	 if(bootstrapValues[j] >= currentThreshold) { 
	 	 		newPredictions[j] <- predictions[j]
	 		}
	 	 else { 
	 			newPredictions[j] <- "unclassified"
	 		}
		}
	 }
	 unclassifiedVector[i] <- length(which(newPredictions == "unclassified"))/length(actual)


# We have the thresholds to filter out based on confidence. 


#	Going through our new predictions.  
	OC <- 0
	MC <- 0

	for(k in seq_along(predictions)) { 
			predicted_Rank <- newPredictions[k]
			actual_Rank <- actual[k]

			if(is.na(predicted_Rank) == FALSE) {

			# When we classify a singleton sequence, were supposed to be getting the classifcationas 
			# "unclassified". If we get anything else, that means the classifier did assignm something to it
			if((actual_Rank %in% singletonGenera)){
				if(predicted_Rank != "unclassified"){

					# This is where we calculate the OCR,
					# So if we come in here, then the sequence we encounter is an 
					# overclassification. So we need to add the predicted rank
					# to a vector.
					OC = OC + 1
					OCpredictions[k] <- predicted_Rank
				}
			}
			# When we try to classifiy a known sequence, we need to take two factors into account. 
			# If the at above that threshold, your thing is given as unclassified, this is because the 
			# confidence for it is lower than that threshold. Futhremore, if it passes the test of having a confidecne
			# that is higher than the threshold we set, then it must pass the test of being classified correctly, which
			# would indicate our misclassification
			else if(predicted_Rank != "unclassified") {
				if(predicted_Rank != actual_Rank){
			
				MC = MC + 1
				MCpredictions[k] <- predicted_Rank
				}
			}
			# cat(actual_Rank %in% knownGenera,"-",predicted_Rank==actual_Rank,"-", MC,'\n')
		}
	}
	# cat('OC = ', OC,'\n')
	OCvector[i] <- OC/length(singletonGenera)
	# cat('OCR =', OCvector[i],'\n')
	MCvector[i] <- MC/lengthOfKnown
	

}

}



sizeOfMCgenera <- table(rdp$genus)[names(table(MCpredictions)[-1])]
sizeOfOCgenera <- table(rdp$genus)[names(table(OCpredictions)[-1])]
sizeofallgenera <- table(rdp$genus)
MCpredictions <- MCpredictions[which(MCpredictions!="")]
OCpredictions <- OCpredictions[which(OCpredictions!="")]
MCpredictions <- unique(MCpredictions)
OCpredictions <- unique(OCpredictions)

# In this case, we can find the location of our location of our genus and find the distances to every
# sequence in the predicted genus.
# Once we do that we can average out the distances and get a number d2.
# Similarly, we can find other sequences in our current genus and find the distances to them from
# our original genus, once we do that we can find the new d1 .
# If we take d1-d2 like this for every misclassification, we can get the histogram with all the distances
