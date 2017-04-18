

#args FILE, CONFIDENCE
args = (commandArgs(TRUE))
library(DECIPHER)
# Load the distance Matrix 
load('/home/adithya/RDP/tfidf/PredictionsRData/allalgorithmsdata/RData/DistanceMatrix.RData')
# Load the file 
load(FILE)
DATASET = as.integer(substr(FILE,1,1))
# Load the trainset
if(DATASET == 1 ) { 
train <- readDNAStringSet('/home/adithya/RDP/tfidf/PredictionsRData/allalgorithmsdata/FastaFiles/trainset16_022016.fa')
} else if (DATASET == 2) { 
train <- readDNAStringSet('/home/adithya/RDP/tfidf/PredictionsRData/allalgorithmsdata/FastaFiles/v4region.fasta')
} else if (DATASET == 3) { 
train <- readDNAStringSet('/home/adithya/RDP/tfidf/PredictionsRData/allalgorithmsdata/FastaFiles/Warcup_v2.fasta')
} 
classes <- sapply(strsplit(names(train), "Root;", fixed=TRUE),
	`[`,
	2L)
# Handle diagonals  
for(i in 1:nrow(d)){ d[i,i] <- 99999}

# Load up the predictions, confidences and the indices
confidence <- eval(parse(text = paste('c',substr(FILE, start =  3, stop = nchar(FILE) - 6), sep = '')))
predicted <- eval(parse(text = paste('p',substr(FILE, start =  3, stop = nchar(FILE) - 6), sep = '')))
index <- eval(parse(text = paste('index',substr(FILE, start =  3, stop = nchar(FILE) - 6), sep = '')))

# Calculate the MC and OC at 80%


t <- table(classes)
weight <- 1/t[classes] # weigh each class equally

singletons <- which(classes %in% names(t[t==1]))
nonsingletons <- which(classes %in% names(t[t > 1]))

# thresholds <- seq(0, 100, 0.1)
# OC <- MC <- CC <- numeric(length(thresholds))

# OC = (# singletons classified) / (# singletons total)
OC <- sum(confidence[singletons] >= CONFIDENCE)/length(singletons)

# MC = (# non-singletons mis-classified) / (# non-singletons classified)
w <- nonsingletons[which(confidence[nonsingletons] >= CONFIDENCE)]
MC <- sum(weight[w[predicted[w] != classes[w]]])/sum(weight[w])

# CC = (# non-singletons classified) / (# non-singletons total)
CC <- sum(weight[w])/sum(weight[nonsingletons])

# w gives us the indices that are misclassified. We will look up those positions
# in the distance matrix d
row_index <- w
col_index <- index[w]   # Grab the classification indices

dis_vector <- numeric(length(row_index))
for(i in seq_along(row_index)){ 
	dis_vector[i] <- d[row_index[i],col_index[i]]	
}

# Save the calculated distance vector.
save(dis_vector, file = paste0(c('MC_dist_',FILE), collapse = ''))
