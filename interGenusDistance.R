# THis script is used  to histogram the distances
# between sequences of the same genus, we get a
# good idea of how far similar sequences are. 


# The way we do this is by taking a submatrix
# form the distance matrix and taking all the distance
# in that submatrix, collapse it and add it to a vector
# or a list and collapse that list and histogram it



# Loading our datafiles
load('distanceMatrix.RData')
load('rdpDataframe.RData')

# Getting the genera

uniqueGenera <- unique(rdp$genus)
genera <- rdp$genus

distList <- list()

for(i in seq_along(uniqueGenera)){
	
	# Get indices
	genusIndices <- which(rdp$genus == uniqueGenera[i])

	# get subset matrixs, this returns a matrix
	subsetMatrix <- d[genusIndices,genusIndices]

	#collapse it into the list
	distList[[i]] <- subsetMatrix[upper.tri(subsetMatrix)]

}

distList <- unlist(distList)

hist(distList)