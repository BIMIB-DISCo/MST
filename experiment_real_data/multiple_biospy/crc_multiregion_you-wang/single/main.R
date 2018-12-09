##############################################################################
###
### MST
###
### Main
###
##############################################################################
### Copyright (c) 2015-2018, The TRONCO Team (www.troncopackage.org)
### email: tronco@disco.unimib.it
### All rights reserved. This program and the accompanying materials
### are made available under the terms of the GNU GPL v3.0
### which accompanies this distribution
##############################################################################

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


## Load dataset

load("dataset_navin_TNBC.RData")


## Perform the inference

source('analysis.R')


## Create onconem

source('generate.onconem.R')

### end of file -- main.R
