library(DECIPHER)
# Function to edit taxid file for leave one out
# Ensure the DECIPHER library is loaded in 
editTaxid <- function(lineNo, option = 'LOOS', trainsetFile, taxidFile) { 
	library(DECIPHER) #Load the appropriate library for parsing
	if(option == 'LOOS') { 

		trainSet <- readDNAStringSet(trainsetFile)
		classes <- sapply(strsplit(names(trainSet), "Root;", fixed=TRUE),`[`,2L)
		table_classes <- table(classes)
		test_seq_class <- classes[lineNo]
		# Only edit the taxid file if it is a singleton
		if(table_classes[test_seq_class] == 1) { 
			
			taxid = parseTaxid(taxidFile) # Parse the taxid file
			classList <- lapply(classes, function(x) { unlist(strsplit(x, ';'))})


			# Find the number of columns to include in our dataframe
			df_col_length <- max(sapply(classList, function(x) { length(x) } )) 

			# Create a dataframe for the training set
			train_df <- data.frame(matrix(0,nrow = length(trainSet),
				 ncol = df_col_length), stringsAsFactors = FALSE)

			for(i in 1:length(trainSet)) {
				g <- classList[[i]] 
				train_df[i, c(1:length(g))] <- g 
			}

			# Now we must lookup our test sequence iteratively
			# to know which lines to discard from our
			# taxid file.

			# get the test sequence
			test_seq <- train_df[lineNo,]
			# trim out any 0s
			test_seq <- test_seq[which(test_seq != 0)]

			# Search every level starting from the lowest to find
			# out how many levels we have to remove from the taxid file

			exitClass <- character(length(test_seq)) #classes that must be removed
			for(i in 1:length(test_seq)) {  # we should traverse from 1 -> length because the genus comes lower down the order than the parent
				class <- unlist(test_seq[i]) # It would help not mix up indices i guess
				if(table(train_df[,i])[class] == 1) {
					exitClass[i] <- class
				}
			}

			# Once all the exit level classes have been obtained, we can now edit
			# the taxid file to remove the levels as displayed in the 

			exitLevel <- which(exitClass != "")
			exitClass <- exitClass[which(exitClass != "")]

			# Now we have level and class index, we can obtain the correct
			# index in the taxid table for the index to remove. 

			exitIndex <- mapply(x = exitClass,y = exitLevel, function(x,y){ 
				return(which(taxid[,'Organismlevelname'] == x & 
						taxid[,'LevelNumber'] == y))					
			})

			# Remove this index from the taxid table and readjust
			# the rest of the sequences
			for(i in 1:length(exitIndex)){ 
					t_exit_lineno <- taxid[exitIndex[i],'LineNumber']
					taxid <- taxid[-exitIndex[i],] # Remove the index in question
					if(length(exitIndex) > 1) {exitIndex[i+1:length(exitIndex)] <- exitIndex[i+1:length(exitIndex)] - 1 }
					# Scale the remaining indices to be removed equally
					taxid[which(taxid[,'LineNumber'] > t_exit_lineno),'LineNumber'] <- taxid[
						which(taxid[,'LineNumber'] > t_exit_lineno),'LineNumber'] - 1	# Decrease the line number
																			# for genera that are greater
																			# than the one being removed. 
					# If the main line is being decreased, this would mean all 
					# parent lines will also have to decrease
					# if they are greater than the one being decreased
					taxid[which(taxid[,'ParentLine'] > t_exit_lineno),'ParentLine'] <- taxid[
						which(taxid[,'ParentLine'] > t_exit_lineno),'ParentLine'] - 1
			}

			# Write a new file with the the remaining data of the taxid
			# file into a new text file that can be read by the 
			# classifier. It must contain the same format as the original
			# taxid file. 

			new_taxid_file <- paste0(c('taxid_',lineNo,'_taxid.txt'), collapse = '')
			sink(new_taxid_file)
			for(i in 1:nrow(taxid)){
				tex <- taxid[i,] 
				tex <- paste0(tex,collapse = '*')
				cat(c(tex,'\n'),append=TRUE)
			 }
			 sink()

		} else { # If it is a non singleton, we will not have to do anything. 
				 #  So we can simply reuse the old taxid file again
			new_taxid_file <- paste0(c('taxid_',lineNo,'_taxid.txt'), collapse = '')
			system(paste('cp',trainsetFile,new_taxid_file,sep=' '))

		}

	}

}


# function to parse the taxid file into a list
parseTaxid <-  function(taxidFile) { 
	taxid = read.table(taxidFile, sep = '*', stringsAsFactors = FALSE)
	colnames(taxid) <- c('LineNumber', 'Organismlevelname', 'ParentLine', 'LevelNumber', 'LevelName')
	return(taxid)
}


# Function to split the training file into train and test. 

train_test_split <- function(trainsetFile, lineNo) { 

	cat("Creating train and test split",'\n')
	train <- readDNAStringSet(trainsetFile) # read the training set. 
	train_t <- train[-lineNo] # split into train 
  	test_t <- train[[lineNo]]	# and test

	data <- names(train) # get the names of the seqeunces	
	trainData <-  data[-lineNo] # get the names of the train data
	testData <- data[lineNo] # get the names of the test data
	trainfile_name <- paste('train_',lineNo,'.fa',sep='')
	testfile_name <- paste('test_',lineNo,'.fa',sep='')
	# You need to remake the train and test files back in fasta
	# format.
	
	## Make the training file 
	sink(trainfile_name)
	for(i in 1:length(trainData)) { 
		cat('>',trainData[i],'\n',sep='')
		cat(as.character(train_t[[i]]),'\n')
	}
	sink()

	## Make the testing file
	sink(testfile_name) 
	cat('>',testData,'\n', sep='')
	cat(as.character(test_t),'\n')
	sink()

}


############### MAIN ###########################
# This is the main command area where we edit the taxid file in order to be able
# to get the correct taxid file. 

# args = (commandArgs(TRUE))


# if(length(args)==0){
#     print("No arguments supplied.")
#     ##supply default values
#     lineNo = 1
#     trainsetFile = 'trainset16_022016.txt'
#     taxidFile = 'trainset16_db_taxid.txt'
# }else{
#     for(i in 1:length(args)){
#          eval(parse(text=args[[i]]))
#     }
# }

# # Execute the editTaxid Function to generate the new taxid file
# editTaxid(lineNo, 'LOOS', trainsetFile, taxidFile)