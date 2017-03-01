

# CalculateOCRMCR takes the filename of the prediction and the confidence
# This should return the a vector with OC and MC for that version.

CalculateOCRMCR <- function(filename_pred, filename_conf){

# filename_pred <- ''
# filename_conf <- ''
load(filename_pred)
load(filename_conf)


load('rdpDataframe.RData')

# In this way, we calculate the MC based on how
# many actually happened in the genus.

# We need to make sure that the actual vector that we compare
# we need to get the indices of the predictions that we made and
# create our bootstrap.


	# Get the indices
	confidencesIndices <- which(confidences != '')
	confidences <- confidences[confidencesIndices]
	predictionsIndices <- which(predictions != '0')
	predictions <- predictions[predictionsIndices]

# First we initialize the data that we need to work with.

thresholdValues <- seq(0.00,1.00,0.01)
predictions <- predictions
confidences <- confidences/100
actual <- rdp$genus[confidencesIndices]
singletonGenera <- names(which(table(actual) == 1))
knownGenera <- names(which(table(actual) != 1))
lengthOfKnown <- length(actual) - length(singletonGenera)
MCvector <- vector(mode = 'double', length = length(thresholdValues))
OCvector <- vector(mode = 'double', length = length(thresholdValues))
bootstrapValues <- confidences
unclassifiedVector <- vector(mode = 'double', length = length(thresholdValues))


# We need to create a dataframe with the following columns
# Actual, Prediction, Unclassified, OC, MC,
# This dataframe has to be set for every new threshold.
# It has 5 columns and length(predictions) rows

classDex <- matrix(0, nrow = length(predictions), ncol = 6)
classDex <- data.frame(classDex)
colnames(classDex) <- c('Actual','Prediction','Confidence','Unclassified','OC','MC')
classDex[,'Actual'] <- actual
# classDex[,'Prediction'] <- predictions
classDex[,'Confidence'] <- confidences

# Start going by each threshold and compare the confidence with teh thresholdvaluess
for(i in seq_along(thresholdValues)) {


	# Calculate the new predictions for based on the theshold we are calculating from
	cat('Threshold = ', thresholdValues[i],'\n')
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

	 # Set our new predictions to the ones that have been modified
	 classDex[,'Prediction'] <- newPredictions
	 # Set our unclassifieds to true or false based on whether that predictions was called
	 # an unclassified or classified.
	 classDex[which(newPredictions == "unclassified"),'Unclassified'] <- TRUE
 	 classDex[which(newPredictions != "unclassified"),'Unclassified'] <- FALSE

 	 # Going through the predictions and comiing up
 	 # with whether that was an MC or OC.

 	 for(k in seq_along(predictions)) {

 	 	 	if(is.na(classDex[,'Prediction'] == FALSE){


	 	 	 	if((classDex[,'Actual'] %in% singletonGenera)){
	 	 	 		if(classDex[,'Unclassified'] == FALSE){
	 	 	 			classDex[,'OC'] <- 1
	 	 	 		}
	 	 	 	}
	 	 	 	else if((classDex[,'Unclassified'] == FALSE)) {
	 	 	 		if(classDex[,'Prediction'] != classDex[,'Actual']){
	 	 	 			classDex[,'MC'] <- 1
	 	 	 		}

	 	 	 	}

			}
 	 }

 	 # The oc vector is calculated as normal.
 	 OCvector[i] <- sum(classDex[,'OC'])/length(singletonGenera)

 	 # Now to calculate the MC's
 	 # In order to do the MC's we need to calculate
 	 # We calculate the MC by taking into account the
 	 # number of MC errors made for that genus present
 	 # from that part of the dataset.

 	 # In this case, we only need to pick out the genera
 	 # Whose MC's are not 0. The genera that are not
 	 # 0 would be the ones where the genus would
 	 # be misclassified

 	 # Scan through the dataframe by genus
 	 # and find out what fraction of the selected
 	 # genus is the genera.

 	 uniqueGenera <- unique(classDex[,'Actual'])
 	 mcv <- sapply(uniqueGenera, FUN = function(x) {
 	 		genSum <- sum(classDex[which(x == classDex[,'Actual']),'MC'])
 	 		genSum <- genSum/length(which(x == classDex[,'Actual']))
 	 		return(genSum)
 	 	})
 	 mcv <- mcv/length(singletonGenera)
 	 MCvector[i] <- mcv
  }
	return(list(mcv = MCvector, ocv = OCvector))
}
