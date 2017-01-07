
args = (commandArgs(TRUE))

if(length(args)==0){
    print("No arguments supplied.")
    ##supply default values
	version=3
	s=32
}else{
    for(i in 1:length(args)){
         eval(parse(text=args[[i]]))
    }
}


plotname <- paste0(c("LOTOv1"), collapse = "")
# jpeg(filename = paste0(c(plotname, ".jpeg"), width = 1366, height = 768))

plot(x = NULL , y = NULL , xlim = c(0.00, 1.00), ylim = c(0.00,1.00), 
    type = 'l' , xlab = '% Classified'
    , ylab = 'Rate', main = 'v3|v3WR|v3WoR|v4WR|v4WoR', xaxs = 'i', yaxs='i', yaxt = 'n')
labels <- c(0, 0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 1.00)
axis(2, labels^0.5,  labels)
axis(2, seq(0,0.1,0.01)^0.5, labels = rep('',11))
axis(2, seq(0.1,1.00, 0.05)^0.5, labels = rep('',19))


# plot(x = NULL , y = NULL , xlim = c(0.00, 1.00), ylim = c(0.00,1.00), 
#     type = 'l' , xlab = 'Threshold'
#     , ylab = 'Rate', main = 'tfidfSintax', xaxs = 'i', yaxs='i')

veks <- c('orig','4mers','5mers','6mers','7mers','8mers','9mers')

for(lab2 in 1:(length(veks) - 1 )) { 

# loadfile1 <- paste(c("../PredictionsRData/",veks[lab2+1],"Confidencesv",version,"s",s,".RData"), collapse = "")
# loadfile2 <- paste(c("../PredictionsRData/",veks[lab2+1],"Predictions.RData"), collapse = "")
# loadfile1 <- 'RDP_V4_region.RData'
loadfile2 <- '8mersPredictions.RData'
loadfile3 <- '8mersConfidencesv3s32.RData'
# load(loadfile1)
load(loadfile3)
load(loadfile2)

# load('origConfidence.RData')
# load('origPredictions.RData')
load('rdpDataframe.RData')



thresholdValues <- seq(0.00,1.00,0.01)


predictions <- predictions
# predictions <- predictions
actual <- rdp$genus
singletonGenera <- names(which(table(rdp$genus) == 1))
knownGenera <- names(which(table(rdp$genus) != 1))
lengthOfKnown <- length(actual)-length(singletonGenera)
MCvector <- vector(mode = 'double', length = length(thresholdValues))
OCvector <- vector(mode = 'double', length = length(thresholdValues))
# bootstrapValues <- confidences
bootstrapValues <- confidences
# if(lab2 == 1){ 
# 			bootstrapValues <- bootstrapValues*100
# 		}
# bootstrapValues <- bootstrapValues/100
unclassifiedVector <- vector(mode = 'double', length = length(thresholdValues))


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
					OC = OC + 1
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

linetype=2
points(x = (1 - unclassifiedVector) , y =  sqrt(OCvector) ,type = 'l', col = 'red' ,lwd = 2, lty = linetype) 

points(x =  (1 - unclassifiedVector) , y =  sqrt(MCvector) ,type = 'l', col = 'blue' ,lwd = 2, lty = linetype) 


}



loadfile1 <- paste(c("../PredictionsRData/",veks[1],"Confidence.RData"), collapse = "")
loadfile2 <- paste(c("../PredictionsRData/",veks[1],"Predictions.RData"), collapse = "")
tryCatch(load(loadfile1), warning = function(x) { 
		cat("File does not exists\n");
		})
tryCatch(load(loadfile2), warning = function(x) { 
		cat("File does not exists\n")
	})
linetype = 1
# load('origConfidence.RData')
# load('origPredictions.RData')
load('rdpDataframe.RData')



thresholdValues <- seq(0.00, 1.00 , 0.01)



predictions <- predictions
# predictions <- predictions
actual <- rdp$genus
singletonGenera <- names(which(table(actual) == 1))
knownGenera <- names(which(table(actual) != 1))
lengthOfKnown <- length(actual)-length(singletonGenera)
MCvector <- vector(mode = 'double', length = length(thresholdValues))
OCvector <- vector(mode = 'double', length = length(thresholdValues))
# bootstrapValues <- confidences
bootstrapValues <- confidences
# if(lab2 == 1){ 
# 			bootstrapValues <- bootstrapValues*100
# 		}
bootstrapValues <- bootstrapValues/100
unclassifiedVector <- vector(mode = 'double', length = length(thresholdValues))


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
					OC = OC + 1
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


points(x = (1 - unclassifiedVector) , y =  sqrt(OCvector) ,type = 'l', col = 'red' ,lwd = 3, lty = linetype) 

points(x =  (1 - unclassifiedVector) , y =  sqrt(MCvector) ,type = 'l', col = 'blue' ,lwd = 3, lty = linetype) 

# points(x = seq(0.00, 1.00, 0.01), y =  unclassifiedVector ,type = 'l', col = 'green' ,lwd = 3, lty = linetype) 


legend("topleft", c('OrigSintax','TfidfSintax4mers','TfidfSintax5mers',
		'TfidfSintax6mers'
		,'TfidfSintax7mers','TfidfSintax8mers', 'TfidfSintax9mers'), lty=c(1,1,2,3,4,5,6), lwd = c(3,2,2,2,2,2,2))



legend(x = locs$x[1],y = locs$y[1], c("OCR","MCR"), 
			 lty = c(1,1), col = c("red","blue"), lwd = c(3,3))
legend(x = locs$x[2], y = locs$y[2], c("Orig","V3"), lty = c(2,3), lwd = c(3,3))

plotloc <- locator()
legend(x = plotloc$x[1], y = plotloc$y[1], c("v5","v3"), lty = c(1,2), lwd = c(2,2))
legend(x = plotloc$x[2], y = plotloc$y[2], c("MCR","OCR"), col = c("blue","red"), lty = c(1,1), lwd = c(2,2))


dev.off()
