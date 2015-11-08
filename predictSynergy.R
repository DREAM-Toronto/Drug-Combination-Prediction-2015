# predictSynergy.R
# Driver for similarity predictions
#
# Purpose: - Makes filenames available in environment
#                so code can be run interactively,
#                from commandline, or source()'d from
#                XvalidatePredictions.R
#          - if arg[5] defines METHOD as "QuickPredict",
#                calls quickSimMat.R to produce DSM and
#                CSM matrices; then calls quickPredict.R
#                to make predictions and write prediction
#                and confidence files for test set.
#          - Reads test data;
#          - Predicts synergy for test data
#          - Writes output
#
# Version: 0.2           
#
# Date:    Nov 7 2015
# Author:  Boris and DREAM team UofT
#          
# V 0.2    Maintenance and refactoring.
# V 0.1    First code
# ==========================================================

setwd(DREAMDIR)


# == CONSTANTS =============================================
#

args <- commandArgs(trailingOnly = TRUE)

# If this script is run interactively - i.e. not source()'d
# or run from commandline, you must uncomment the
# assignments below and give them the filenames you need.
# args <- NULL
# args[1] <- "../validate/ch1_train_combination_and_monoTherapy.csv"
# args[2] <- "../validate/ch1_leaderBoard_monoTherapy.csv"      
# args[3] <- "../validate/08.01-prediction.csv"
# args[4] <- "../validate/08.01-combination_priority.csv"
# args[5] <- "QuickPredict"


TRAINING_SET_FILE <- args[1]
TEST_SET_FILE     <- args[2]
PREDICTED_FILE    <- args[3]
PRIORITY_FILE     <- args[4]
METHOD            <- args[5]   

# == PACKAGES ==============================================


# == FUNCTIONS =============================================
#


# == MAIN ==================================================
#

if (METHOD == "QuickPredict") {
	
# >>> COMPUTE SIMILARITIES
# >>> source the code that computes
# >>> compound and cell similarities and
# >>> creates DSM and CSM matrices
source("quickSimMat.R")

# >>> MAKE PREDICTIONS
# >>> source the code that makes the actual
# >>> predictions.
source("quickPredict.R")

# Done

}



# [END]






