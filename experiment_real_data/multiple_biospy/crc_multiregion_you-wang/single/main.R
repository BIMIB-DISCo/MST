if (!require(TRONCO)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("TRONCO")
    library(TRONCO)
}

if (!require(oncoNEM)) {
	library(devtools)
	install_bitbucket("edith_ross/oncoNEM")
	library(oncoNEM)
}

options(scipen = 999)

#load dataset
load("dataset_navin_TNBC.RData")

# perform the inference
source('analysis.R')

# create onconem
source('generate.onconem.R')