#!/bin/bash
#THis script calls the plotter with the version of the code and the bootstrap size 
# They must exist in the PredictedRData Folder.
Rscript plot.R version=$1 s=$2
