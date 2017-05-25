

# Function to edit the training data
# train : training data to create new modified training file from
# retentionPercent : what percent of the training data should be retained
# addRandomFunny : Adds random funny sequences inside just to mess with the classifier. 
# modificationPercent :  What percent of the seqeunces would you like to modify at maximum
activateMutation <- function(train, retentionPercent, mutationPercent) { 

		# Get the indices that we want to modify 
		retentionIndices <- sample(length(train), round(retentionPercent*length(train)))

		# get the subset at the indices
		modtrain <- train[retentionIndices]
		# Add insertions of random between A,T,G,C
		modifiedSequences <- sapply(as.character(modtrain), function(x) { 
				
				# option <- sample(5,1)

				# if(option == 1) { 
				# 		# Indels
				# } else if (option == 2) { 
				# 		# Deletions
				# } else if (option == 3) { 
				# 		# Reverse complement a subsequence
				# } else if (option == 4) { 
				# 		# Swap two subsets 
				# }

				# Pick a random position to split
				splitPos <- sample(nchar(x),1)

				# Pick the nucleotide to clone at this position
				splitDNA <- unname(unlist(strsplit(x, '')))
				nucleotide <- splitDNA[splitPos]

				# Clone this nucleotide a buncha times 
				# based on mutationPercent
				mutation <- unname(rep(nucleotide, nchar(x)*mutationPercent))

				# Insert this mutation at the position of the  split
				mutatedDNA <- paste0(c(splitDNA[1:splitPos-1], mutation, splitDNA[splitPos:length(splitDNA)][-1]), collapse = '')

				# Return the mutatedDNA
				return(mutatedDNA)
			})

		# Create the file to store the new seqeunces
		FILE = 'mutatedTestData.fa'
		sink(FILE) 
		# store the new values inside the new file 
		for(i in seq_along(modifiedSequences)) { 
			cat(paste0(c('>TestSequence', i), collapse = ''),'\n')
			cat(modifiedSequences[[i]],'\n')
		}
		# close 
		sink()
		return(FILE)
}