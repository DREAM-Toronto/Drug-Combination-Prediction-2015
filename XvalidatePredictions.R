# XvalidatePredictions.R
#
# Purpose: From input Training Data, select a subset for
#              Training and prediction, run
#              similarity computations,
#              predict synergy scores, and evaluate
#              success.
#            
#
# Date:    Nov 5 2015
# Author:  Boris and DREAM team UofT
#          
# V 0.1    First code
# ==========================================================

setwd(DREAMDIR)


# == CONSTANTS =============================================
#
FRAC_HOLD <- 0.1  # Fraction of holdout data

#
EXCLUDE_POOR_QA <- FALSE  


# 
INFILE <- "../Challenge Data/Drug Synergy Data/ch1_train_combination_and_monoTherapy.csv"

OUTFILE_TRAIN <- "../validate/xValTrain.csv"

OUTFILE_HOLD <- "../validate/xValTest.csv"

TEST_PREDICTED_FILE <- "../validate/test_predicted.csv"


# == PACKAGES ==============================================


# == FUNCTIONS =============================================
#


# == MAIN ==================================================
#

# base file from which to select
base <- read.csv(INFILE, stringsAsFactors=FALSE)

if (EXCLUDE_POOR_QA) {
	base <- base[base[,"QA"] == 1, ]
}

# split "base" into "training" and "holdout" data 
nHold <- round(nrow(base) * FRAC_HOLD) # number of holdouts
iHold <- sample(1:nrow(base), nHold)   # row index of holdouts

trainSet <- base[-(iHold), ]
holdSet  <- base[iHold, ]
holdSet[, SYNERGY_SCORE] <- NA  # remove synergy values from
                                # holdout data

write.csv(trainSet, OUTFILE_TRAIN)
write.csv(holdSet, OUTFILE_HOLD)

# == RUN PREDICTION =======

TRAIN_FILE <- OUTFILE_TRAIN
TEST_FILE <- OUTFILE_HOLD

# predictSynergy.R is written to take input training and test
# filenames from commandline. If it it source()'ed instead
# redefine the commandArgs() function as below. 

commandArgs <- function(trailingOnly) {
	return(c(TRAIN_FILE, TEST_FILE, TEST_PREDICTED_FILE))
}
source("predictSynergy.R")

# predictSynergy.R writes its result into 
# TEST_PREDICTED_FILE

prediction <- read.csv(TEST_PREDICTED_FILE, stringsAsFactors=FALSE)

predSynScores <- prediction[, SYNERGY_SCORE]
trueSynScores <- base[iHold, SYNERGY_SCORE]

plot(trueSynScores, predSynScores)

# calculate correlation
cor(trueSynScores, predSynScores)




# [END]
