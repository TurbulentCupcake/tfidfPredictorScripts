
args = (commandArgs(TRUE))

if(length(args)==0){
    print("No arguments supplied.")
    ##supply default values
}else{
    for(i in 1:length(args)){
         eval(parse(text=args[[i]]))
    }
}

if(test == TRUE) { 
 	cat('True')
 } else  {
	cat('Fals')
}