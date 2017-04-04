# Definitive script to load and calculate and save a vector
# that contains the unbalanced/balanced ratio 

library(DECIPHER)
TRAINFILE <- "trainset16_022016.fa"
ISCOMBINATION <- TRUE

if(ISCOMBINATION) {
	# Specify your combination of break in slope and shift
	slope <- c(1,100,10)
	shift <- c(1,100,1)

	for(i in seq_along(slope)){ 
		for(j in seq_along(shift)){ 
			loadfilebalanced <- paste0(c('balanced_',slope,'_',shift,'.RData'), collapse = '')
			loadfileunbalanced <- paste0(c('unbalanced_',slope,'_',shift,'.RData'), collapse = '')
			ratio <- calculateRatioBalUnbal(loadfileunbalanced, loadfilebalanced)
			savefile <- paste0(c('ratio_',slope,'_',shift,'.RData'), collapse = '')

		}
	}
} else { 
	# specify the name of the load file
	loadfilebalanced <- 'tophitBalanced.RData'
	loadfileunbalanced <- 'tophitUnbalanced.RData'
	ratio <- calculateRatioBalUnbal(loadfilebalanced, loadfileunbalanced)
}

calculateRatioBalUnbal <- function(loadfileunbalanced, loadfileunbalanced){ 
 

}


