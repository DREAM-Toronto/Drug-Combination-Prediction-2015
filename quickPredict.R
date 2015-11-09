# quickPredict.R
# Quick calculation of predictions for target data.
#
# Calculates a synergy score for each test set 
# drug/drug/cell combination, from an average of observed
# synergy scores in the training set, weighted by the
# similarity of the resp. drugs and cells in the tranining
# set to the test set combination. 
#
# Version: 0.2.1           
#
# Date:    Nov 8 2015
# Author:  Boris and DREAM team UofT
#          
# V 0.2.1  Update handling of Bad Rows.
#          Add progress messages.
# V 0.2    Maintenance and refactoring.
#          Bugfix in calculating CSM
# V 0.1    First code
# ==========================================================

# This file needs to be source()'d from predictSynergy.R 


# TEST_FILE <- "../validate/xValTest.csv"
# PREDICT_FILE <- "../validate/xValTest_predicted.csv"


# == CONSTANTS =============================================
#


TOP_N <- 5  # Use top N best similarities to calculate
            # weighted average for similarity score


# == FUNCTIONS =============================================
#

checkABC <- function(dat, drugs, cells) {
    # Returns index of rows in a set that contain drugs or
    # cells that are not in the list of drugs and cells that 
    # were used to compile the DSM and CSM.
	# "drugs" contains drug names
	# "cells" contains cell names
    
    x <- NULL
    x <- c(x, which(!(dat[ , "COMPOUND_A"] %in% drugs)))
    x <- c(x, which(!(dat[ , "COMPOUND_B"] %in% drugs)))
    x <- c(x, which(!(dat[ , "CELL_LINE"]  %in% cells)))

    return(unique(x))
}



# == MAIN ==================================================
#

if (VERBOSE) {cat(paste("    Predicting synergy ..."))}

testData <- read.csv(TEST_SET_FILE,
                     header = TRUE,
                     row.names = NULL,
                     stringsAsFactors = FALSE)

# If DSM and CSM have not just been compiled, we need to 
# read them from somewhere now ...
# DSM <- read.csv(DSM_FILE, stringsAsFactors=FALSE)
# CSM <- read.csv(CSM_FILE, stringsAsFactors=FALSE)

badRows <- checkABC(testData, drugs, cells)
if (length(badRows) > 0) {
	# during xValidation testing, it could happen that
	# the holdout set contains drugs or cells that were
	# not present in the training set. This should not
	# happen in the real submissions, and it should be
	# infrequent enough that we can simply remove the
	# offending row(s) from the test data.
	
	testData <- testData[-(badRows),]
	if (VERBOSE) {
		cat(" removed Bad Row(s) ")
		cat(paste(badRows, collapse = " "))
		cat(" ...")
	}

}
nTest <- nrow(testData)


# Initialize simSyn to store holds the probability of 
# similarity to arequested unknown combination and the
# corresponding  similarity score, for each combination 
# in the training data.
simSyn <- matrix(0, nrow=nTrain, ncol=2)
colnames(simSyn) <- c("pSim", "SYN")

# == PREDICT SIMILARITY SCORES =======

for (iTest in 1:nTest) {
	A_test <- testData[iTest, "COMPOUND_A"]
	B_test <- testData[iTest, "COMPOUND_B"]
	C_test <- testData[iTest, "CELL_LINE"]
	for (iTrain in 1:nTrain) {
        dA <- DSM[A_test, trainingData[iTrain, "COMPOUND_A"]]
        dB <- DSM[B_test, trainingData[iTrain, "COMPOUND_B"]]
        dC <- CSM[C_test, trainingData[iTrain, "CELL_LINE"]]
        simSyn[iTrain, "pSim"] <- dA * dB * dC
        simSyn[iTrain, "SYN"] <- trainingData[iTrain, "SYNERGY_SCORE"]
	}
	# order descending by pSim
	simSyn <- simSyn[order(simSyn[,"pSim"], decreasing=TRUE), ]
	
	# calculate weighted average for TOP_N rows
	avSYN <- simSyn[1:TOP_N, "pSim"] * simSyn[1:TOP_N, "SYN"]
	sumSYN <- sum(simSyn[1:TOP_N, "pSim"])
	if (sumSYN == 0) { # all probabilities zero
		avSYN <- 0
	} else {
        avSYN <- sum(avSYN) / sum(simSyn[1:TOP_N, "pSim"])
	}
	
	# store result
	testData[iTest, "SYNERGY_SCORE"] <- avSYN
}


# check if all synergy scores have been predicted
if (any(is.na(testData[,"SYNERGY_SCORE"]))) {
	stop("PANIC: remaining NAs in testData.")
}

# rename "SYNERGY_SCORE" to "PREDICTION"
colnames(testData) <- sub("SYNERGY_SCORE", "PREDICTION", colnames(testData))

# == WRITE PREDICTION TO FILE =======
# write output in challenge 1 "submission format
write.csv(testData[ , c("CELL_LINE", "COMBINATION_ID", "PREDICTION")],
          PREDICTED_FILE,
          row.names = FALSE)


# == CALCULATE QUICK CONFIDENCE ===========
# "confidence" is the belief that a drug combination will be
# synergistic. As a quick approach, we calculate confidence
# as the scaled average rank of a combination's synergy
# score: the higher the predicted synergy score, the stronger
# our belief that there is synergy.

conf <- data.frame("COMBINATION_ID" = sort(unique(testData[,"COMBINATION_ID"])), 
                   "CONFIDENCE" = 0,
                   stringsAsFactors=FALSE)
# head(confidence)

for (i in 1:nrow(conf)) {
	conf[i,"CONFIDENCE"] <- mean(testData[testData[ ,"COMBINATION_ID"] == conf[i,"COMBINATION_ID"], "PREDICTION"])
}
# head(conf)

conf[ ,"CONFIDENCE"] <- rank(conf[,"CONFIDENCE"])/nrow(conf)
# head(conf)

# == WRITE CONFIDENCE TO FILE =======
write.csv(conf,
          PRIORITY_FILE,
          row.names = FALSE)

# Done.

if (VERBOSE) {cat(paste(" Done\n"))}


# [END]
