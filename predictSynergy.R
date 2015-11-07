# predictSynergy.R
#
# Purpose: - Reads Training Data
#          - Gets compound and cell line  similarity;
#          - Reads test data;
#          - Predicts synergy for test data
#          - writes output
#            
#
# Date:    Nov 6 2015
# Author:  Boris and DREAM team UofT
#          
# V 0.1    First code
# ==========================================================

setwd(DREAMDIR)


# == CONSTANTS =============================================
#

args <- commandArgs(trailingOnly = TRUE)

# If this script is run interactively - i.e. not source()'d
# or run from commandline, you must uncomment the
# assignments below and give them the filenames you need.

# args <- character(3)
# args[1] <- "../validate/xValTrain.csv"
# args[2] <- "../validate/xValTest.csv"      
# args[3] <- "../validate/test_predicted.csv"

TRAIN_FILE <- args[1]
TEST_FILE  <- args[2]
PRED_FILE  <- args[3]


# == PACKAGES ==============================================


# == FUNCTIONS =============================================
#


# == MAIN ==================================================
#

trainingData <- read.csv(TRAIN_FILE, stringsAsFactors=FALSE)


# >>> COMPUTE SIMILARITIES
# >>> here we need to source the code that computes
# >>> compound and cell similarities


# >>> Let's store the results in two files that we will
# >>> assign to DSM_FILE (Drug Similarity Matrix) and
# >>> CSM_FILE (Cell Similarity Matrix)


# >>> MAKE PREDICTIONS
# >>> Here we need to source the code that makes the actual
# >>> predictions.

system(paste("makePrediction.py",
             "TEST_FILE",
             "DSM",
             "CSM",
             "PRED_FILE"))

# Done: PRED_FILE has been written


# [END]






