


	# source('CalculateOCRMCR.R')

# # Step 1, calculate the OCR MCR by going through each vector. 
# confidenceNames <- c('Sintax_Confidence_k16_lineLength100',
# 'Sintax_Confidence_k16_lineLength1000',
# 'Sintax_Confidence_k16_lineLength250',
# 'Sintax_Confidence_k16_lineLength50',
# 'Sintax_Confidence_k16_lineLength500',
# 'Sintax_Confidence_k24_lineLength100',
# 'Sintax_Confidence_k24_lineLength1000',
# 'Sintax_Confidence_k24_lineLength250',
# 'Sintax_Confidence_k24_lineLength50',
# 'Sintax_Confidence_k24_lineLength500',
# 'Sintax_Confidence_k32_lineLength100',
# 'Sintax_Confidence_k32_lineLength1000',
# 'Sintax_Confidence_k32_lineLength250',
# 'Sintax_Confidence_k32_lineLength50',
# 'Sintax_Confidence_k32_lineLength500',
# 'Sintax_Confidence_k40_lineLength100',
# 'Sintax_Confidence_k40_lineLength1000',
# 'Sintax_Confidence_k40_lineLength250',
# 'Sintax_Confidence_k40_lineLength50',
# 'Sintax_Confidence_k40_lineLength500',
# 'Sintax_Confidence_k8_lineLength100',
# 'Sintax_Confidence_k8_lineLength1000',
# 'Sintax_Confidence_k8_lineLength250',
# 'Sintax_Confidence_k8_lineLength50',
# 'Sintax_Confidence_k8_lineLength500')
	
# predictionNames <- c('Sintax_Predictions_k16_lineLength100',
# 'Sintax_Predictions_k16_lineLength1000',
# 'Sintax_Predictions_k16_lineLength250',
# 'Sintax_Predictions_k16_lineLength50',
# 'Sintax_Predictions_k16_lineLength500',
# 'Sintax_Predictions_k24_lineLength100',
# 'Sintax_Predictions_k24_lineLength1000',
# 'Sintax_Predictions_k24_lineLength250',
# 'Sintax_Predictions_k24_lineLength50',
# 'Sintax_Predictions_k24_lineLength500',
# 'Sintax_Predictions_k32_lineLength100',
# 'Sintax_Predictions_k32_lineLength1000',
# 'Sintax_Predictions_k32_lineLength250',
# 'Sintax_Predictions_k32_lineLength50',
# 'Sintax_Predictions_k32_lineLength500',
# 'Sintax_Predictions_k40_lineLength100',
# 'Sintax_Predictions_k40_lineLength1000',
# 'Sintax_Predictions_k40_lineLength250',
# 'Sintax_Predictions_k40_lineLength50',
# 'Sintax_Predictions_k40_lineLength500',
# 'Sintax_Predictions_k8_lineLength100',
# 'Sintax_Predictions_k8_lineLength1000',
# 'Sintax_Predictions_k8_lineLength250',
# 'Sintax_Predictions_k8_lineLength50',
# 'Sintax_Predictions_k8_lineLength500')

# This function returns the MC and OC vector for a given percentage across all the data

rangeofs <- seq(1,16,1)
predictionNames <- sapply(rangeofs, FUN = function(x) { 
	return(paste(c(x,'_50_predictionv2.RData'), collapse = ''))
})
confidenceNames <- sapply(rangeofs, FUN = function(x){  
	return(paste(c(x,'_50_confidencev2.RData'), collapse = ''))
})

isopercent <- function(percent = 0.80, predictionNames, confidenceNames, indexFile){ 


	# create the memory to hold the data

	cat('Isolating percentage value = ', percent,'\n')
	ocatpercent <- vector(mode = 'double', length = length(predictionNames))
	mcatpercent <- vector(mode = 'double', length = length(predictionNames))

	cat('Begin Calculation Process .. ','\n')
	# calculate the OCR and MCR for each data vector
	for(i in seq_along(predictionNames)){
		cat('Calculating OCR and MCR for ',predictionNames[i],'\n')
		filename_pred <- predictionNames[i]
		filename_conf <- confidenceNames[i]
		# Get the vector of OCs and MCs
		OCMCUnc <- CalculateOCRMCR(filename_pred, filename_conf, indexFile)
		# From both the vectors, get the percentage that we want and store it in a vector. 
		index <- as.integer(percent*100)
		cat('Getting the oc and mc rate at ',percent,'\n')
		classified <- 1 - OCMCUnc$uncl
		# in order to extrapolate, we need to find out if our percent can actually
		# exist as a percent classified within the classified vector calculated above.
		if ((percent < min(classified)) || (percent > max(classified))) {
			ocatpercent[i] <- NA
			mcatpercent[i] <- NA 

		} else  { 
			# If it isnt less than the minimum, it means we can extrapolate from datapoints
			# given classified vector. 
			index <- which.min(abs(classified-percent))
			# Now we have the percent index th
			yaxisOC <- c(OCMCUnc$ocv[index], OCMCUnc$ocv[index-1])
			yaxisMC <- c(OCMCUnc$mcv[index], OCMCUnc$mcv[index-1])
			xaxis <- c(classified[index], classified[index-1])

			# We calculate the percent classified first for OC
			# in order to do this, we find the slope first, which
			# is given by rise over run. y=mx + 
			slope_OC <-  (yaxisOC[2] - yaxisOC[1])/(xaxis[2] - xaxis[1])
			b_OC <- yaxisOC[2] - slope_OC*xaxis[2]
			ocatpercent[i] <- slope_OC*percent + b_OC

			# Do the same thing for MC
			slope_MC <- (yaxisMC[2] - yaxisMC[1])/(xaxis[2] - xaxis[1])
			b_MC <- yaxisMC[2] - slope_MC*xaxis[2]
			mcatpercent[i] <- slope_MC*percent + b_MC
			

		}


		cat('OC at ----- ',percent,'% = ',ocatpercent[i],'\n')
		cat('MC at ----- ',percent,'% = ',mcatpercent[i],'\n')
	}


		nameVec <- lapply(predictionNames, FUN = function(x){
			subpart <- substr(x, 20, nchar(x))
			t1 <- strsplit(subpart, '_')
			# t1[1] <- strsplit(t1[1],'')[2]
			# t1[2] <- strsplit(t1[2],'')[2]
			# t1 <- paste(c(t1[1],'_',t1[2]),collapse = '')
			return(t1[[1]])
		})
		nameVec2 <- lapply(nameVec, FUN = function(x) { 
			s1 <- unlist(strsplit(x[1], 'k'))[2]
			s2 <- unlist(strsplit(x[2], 'lineLength'))[2]
			finname <- paste(c(s1,'_',s2), collapse = '')
			return(finname)
		})


	# Now that we recieve the vectors with the values at a given percentage, we have to index
	# them correctly in order to be able to identify and plot them correc

	names(mcatpercent) <- nameVec2
	names(ocatpercent) <- nameVec2
	
	cat('Returning vectors','\n')
	return(list(mcvp = mcatpercent, ocvp = ocatpercent))

}

# This function get the MCs and OCs for plotting by prepping data
# to be represented by length and by subsample(s) value
getSnLnPlot <- function(mcvocvlist){ 


	# In this function, we get the vectors from the isolatepercent
	# mcvocvlist contains the mcatpercent and ocatpercent vectors. 

	# Parse the names from the vectors such that we get the vectors in terms of 
	# the subsample values and length value
	mcvp <- mcvocvlist$mcvp
	ocvp <- mcvocvlist$ocvp
	# Get the names
	id <- names(mcvp)
	# string split all the names 
	splt <- sapply(id, FUN = function(x) { strsplit(x,'_')})

	# Now that we have the splits, we need to store them in a matrix based on their sizes

	colsID <- unique(sapply(splt, FUN = function(x){x[[1]]}))
	rowsID <- unique(sapply(splt, FUN = function(x){x[[2]]}))

	mcvpMatrix <- matrix(mcvp, nrow = length(rowsID), ncol = length(colsID), dimnames = list(rowsID, colsID))
	ocvpMatrix <- matrix(ocvp, nrow = length(rowsID), ncol = length(colsID), dimnames = list(rowsID, colsID))	

	plot(x = NULL, y = NULL, ylim = c(0.00,1.00), xlim = c(1,length(rowsID)), xaxt = 'n',
		xlab = 'Sequence Lengths', ylab = 'Rate', yaxt = 'n')
	labels <- c(0, 0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 1.00)
	axis(2, labels^0.5,  labels)
	axis(2, seq(0,0.1,0.01)^0.5, labels = rep('',11))
	axis(2, seq(0.1,1.00, 0.05)^0.5, labels = rep('',19))
	axis(1, c(1:length(rowsID)), labels = rowsID)

	for(i in 1:length(colsID)) { 
		# plot the misclassification
		points(x = c(1:length(rowsID)), y = sqrt(mcvpMatrix[,colsID[i]]), type = 'l', lty = i, lwd = 2, col = 'blue')
		# plot the overclassification
		points(x = c(1:length(rowsID)), y = sqrt(ocvpMatrix[,colsID[i]]), type = 'l', lty = i, lwd = 2, col = 'red')
	}


	plotloc <- locator()
	legend(x = plotloc$x[1], y = plotloc$y[1], c(colsID), lty = c(1:length(colsID)), lwd = c(rep(2,length(colsID))))
	legend(x = plotloc$x[2], y = plotloc$y[2], c("MCR","OCR"), col = c("blue","red"), lty = c(1,1), lwd = c(2,2))



	plot(x = NULL, y = NULL, ylim = c(0.00,1.00), xlim = c(1,length(rowsID)), xaxt = 'n',
		xlab = 'Replicate Size', ylab = 'Rate', yaxt = 'n')
	labels <- c(0, 0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 1.00)
	axis(2, labels^0.5,  labels)
	axis(2, seq(0,0.1,0.01)^0.5, labels = rep('',11))
	axis(2, seq(0.1,1.00, 0.05)^0.5, labels = rep('',19))
	axis(1, c(1:length(colsID)), labels = colsID)

	for(i in 1:length(rowsID)) { 
		# plot the misclassification
		lines(x = c(1:length(colsID)), y = sqrt(mcvpMatrix[rowsID[i],]), type = 'l',lty = i, lwd = 2, col = 'blue')
		# plot the overclassification
		lines(x = c(1:length(colsID)), y = sqrt(ocvpMatrix[rowsID[i],]), lty = i, lwd = 2, col = 'red')
	}


	plotloc <- locator()
	legend(x = plotloc$x[1], y = plotloc$y[1], c(rowsID), lty = c(1:length(rowsID)), lwd = c(rep(2,length(rowsID))))
	legend(x = plotloc$x[2], y = plotloc$y[2], c("MCR","OCR"), col = c("blue","red"), lty = c(1,1), lwd = c(2,2))



}
