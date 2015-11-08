# XvalidatePredictions.R
#
# Purpose: From input Training Data, select a subset for
#              Training and prediction, run
#              similarity computations,
#              predict synergy scores, and evaluate
#              success.
#
# Version: 0.2           
#
# Date:    Nov 7 2015
# Author:  Boris and DREAM team UofT
#          
# V 0.2    Maintenance and refactoring.
#          Loop over Xval runs.
#          Calculate random correlations
# V 0.1    First code
# ==========================================================

setwd(DREAMDIR)
source("DREAMutilities.R")

# == CONSTANTS =============================================
#
METHOD <- "QuickPredict"

RUN_PATH <- "../validate/" # Directory for validation run data
RUN_ID <- "07.01" # Filename prefix

FRAC_HOLD <- 0.1  # Fraction of holdout data

N_RUNS <- 1       # Number of cross-validation runs

EXCLUDE_POOR_QA <- TRUE  # if TRUE, ignore all experiments
                         # with QA != 1

# Probably no need to change below here ====================                         
# Master file that contains monotherapy data
MASTER_DATA_FILE <- "../Challenge Data/Drug Synergy Data/ch1_train_combination_and_monoTherapy.csv"

# Training subset filename
TRAINING_SET_FILE <- sprintf("%s%s_TrainingSet.csv", RUN_PATH, RUN_ID)

# Holdout subset filename
TEST_SET_FILE <- sprintf("%s%s_HoldoutSet.csv", RUN_PATH, RUN_ID)

# Prediction results
PREDICTED_FILE <- sprintf("%s%s_prediction.csv", RUN_PATH, RUN_ID)
PRIORITY_FILE <- sprintf("%s%s_combination_priority.csv", RUN_PATH, RUN_ID)


# == PACKAGES ==============================================


# == FUNCTIONS =============================================
#


# == MAIN ==================================================
#

# master file from which to select
master <- read.csv(MASTER_DATA_FILE, stringsAsFactors=FALSE)

if (EXCLUDE_POOR_QA) {
	master <- master[master[,"QA"] == 1, ]
}

allCorrelations <- NULL
for (iRun in 1:N_RUNS) {
	
	# == MAKE TRAINING AND HOLDOUT SETS =======
	nHold <- round(nrow(master) * FRAC_HOLD) # number of holdouts
	iHold <- sample(1:nrow(master), nHold)   # random row index of holdouts
	
	trainSet <- master[-(iHold), ]
	holdSet  <- master[iHold, ]
	holdSet[, "SYNERGY_SCORE"] <- NA  # remove synergy values from
	                                  # holdout data
	
	write.csv(trainSet, TRAINING_SET_FILE, row.names=FALSE)
	write.csv(holdSet, TEST_SET_FILE, row.names=FALSE)
	
	# == RUN PREDICTION =======
	
	# predictSynergy.R is written to take input training and test
	# filenames from commandline. If it it source()'ed instead
	# redefine the commandArgs() function as below. 
	
	commandArgs <- function(trailingOnly) {
		return(c(TRAINING_SET_FILE,
		         TEST_SET_FILE,
		         PREDICTED_FILE,
		         PRIORITY_FILE,
		         METHOD))
	}
	source("predictSynergy.R")
	
	# predictSynergy.R writes its result into 
	# PREDICTED_FILE and PRIORITY_FILE
	
	# == READ AND COMBINE RESULTS =======
	pred <- read.csv(PREDICTED_FILE, stringsAsFactors=FALSE)
	conf <- read.csv(PRIORITY_FILE, stringsAsFactors=FALSE)
    # head(pred)
    # head(conf)
	
	results <- data.frame(RUN = rep(iRUN, nrow(pred)),
	                      CELLS  = pred[ , "CELL_LINE"],
	                      ID  = pred[ , "COMBINATION_ID"],
	                      TRUE_SYN  = master[iHold, "SYNERGY_SCORE"],
	                      PRED_SYN  = pred[ , "PREDICTION"],
	                      SYN_CONF  = rep(0, nrow(pred)),
	                      stringsAsFactors = FALSE)
    for (i in 1:nrow(results)) { # add confidence values to results
    	results[i, "SYN_CONF"] <- conf[ (conf[,"COMBINATION_ID"] == results[i, "ID"]) ,
    	                                "CONFIDENCE" ]
    }
    # head(results)
    
    if (iRun == 1) {
    	allResults <- results
    } else {
    	allResults <- rbind(allResults, results)
    }
    allCorrelations[iRun] <- cor(results[,"TRUE_SYN"],
                                 results[,"PRED_SYN"],
                                 use="complete.obs")
}   # end for (iRun in 1:N_RUNS)


# plot true scores vs. predictions
true <- allResults[,"TRUE_SYN"]
pred <- allResults[,"PRED_SYN"]
plot(true, pred,
     xlim = c(min(c(true, pred)), max(c(true, pred))),
     ylim = c(min(c(true, pred)), max(c(true, pred))),
     xlab = "True synergy scores",
     ylab = "Predicted synergy scores",
     main = sprintf("%s crossvalidation with %d%% holdouts, %d runs",
                    METHOD,
                    FRAC_HOLD * 100,
                    N_RUNS),
     cex.main = 0.8
    )
abline(h = 0, lwd = 0.5, col = "#CCCCCC")
abline(v = 0, lwd = 0.5, col = "#CCCCCC")


# get distribution of random correlations
nRandRuns <- 10000  # Nice. Takes about 1 second...
corRand <- numeric(nRandRuns)

for (i in 1:nRandRuns) {
	corRand[i] <- cor(true, sample(pred, length(pred)))
}

hist(corRand, 
     breaks=20, 
     xlim = c(min(min(corRand), min(allCorrelations)) * 1.2,
              max(max(corRand), max(allCorrelations)) * 1.2)
    )
for (i in 1:length(allCorrelations)) {
	abline(v = allCorrelations[i], col="#AA0000")
}


# [END]
