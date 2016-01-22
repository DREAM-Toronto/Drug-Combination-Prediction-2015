# DREAM_Combi_main.R
#
# Purpose: Configure and drive a prediction run for the Drug
#             Combination Prediction Challenge. This is a 
#             framework to tie together our code assets.
#
# Version: 0.1           
#
# Date:    Jan 22 2016
# Author:  Boris and DREAM team UofT
#
# V 0.1    First code
#
# TODO:    
#
# ==========================================================

# == DOCUMENTATION =========================================

# == PURPOSE ==
# Note: Provide a summary of the purpose of a particular
#        run here, then save a version to keep a record.

# == RESULTS ==
# Note: Record key results here.


# == SETUP =================================================

setwd(DREAMDIR)
DCP <- new.env() # Drug Combination Prediction "environment"

# == CONFIGURE =============================================

DCP$PATH <- "../test.01/" # path for run output
DCP$PREFIX <- "01.test"  # filename prefix for this run

# == HELPERS ==
source("DCP_HelperFunctions.R") 

# == INITIALIZE ============================================

# Enter the correct source files for your actual 
# functions here:

ABCSLOADER <- "loadABCSFromMaster.0.1.R"            # functions to load synergy score tables
FEATUREMAKER <- "makeFeaturesFromMonotherapy.0.1.R" # functions to make feature tables
FEATURECOMBINER <- "combineDC_Features.0.1.R"       # functions to combine features for training
PREDICTOR <- "quickPredictFunctions.0.1.R"          # functions to make the actual prediction
FORMATTER <- "formatSubmissionFunctions.0.1.R"      # functions to reformat a prediction into a submission
SCORER <- "scoreSubmissionFunctions.0.1.R"          # functions to evaluate prediction scores


# == I: FEATURES ===========================================


# == DEFINE TRAINING DATA AND TEST DATA ==

# No precondition.

# EITHER make a random training and test set (Mode: RANDOM) ...
DCP$MODE <- "RANDOM"
DCP$HOLDOUT_FRACTION <- 0.33

# OR: read ABCS training and test sets from file (Mode: FILE) ...
# DCP$MODE <- "FILE"
# DCP$ABCS_TRAINING_FILE <- "myAwesomeTrainingSet.csv"
# DCP$ABCS_TEST_FILE     <- "myAwesomeTestSet.csv"

# OR: read all available training data and the
#     challenge submission test set (Mode: CHALLENGE).
# DCP$MODE <- "CHALLENGE"
# DCP$HOLDOUT_FRACTION <- 0
# DCP$ABCS_TEST_FILE     <- "leaderboard.csv"

if (DCP$MODE == "RANDOM") {

    source(ABCSLOADER)

    # Load A-B-C-Synergies
    # The return value of loadABCSmaster() is a dataframe
    # with four columns: 
    #    Drug A name
    #    Drug B name
    #    Cell line name
    #    Synergy score

    # set.seed(112358)   # uncomment for reproducible runs
                         # for testing/debugging purposes only!
    ABCS_Master <- loadABCSmaster()
    nHold <- round(nrow(ABCS_Master) * DCP$HOLDOUT_FRACTION) # number of holdouts
    iHold <- sample(1:nrow(ABCS_Master), nHold)   # random row index of holdouts
    
    ABCS_Training <- ABCS_Master[-(iHold), ]
    ABCS_Test  <- ABCS_Master[iHold, ]

    # Uncomment if you want to store the results
    # TrainOutFile <- sprintf("%s%s_TrainingSet.csv", DCP$PATH, DCP$PREFIX)
    # TestOutFile  <- sprintf("%s%s_TestSet.csv", DCP$PATH, DCP$PREFIX)
    # write.csv(ABCS_Training, TrainOutFile, row.names=FALSE)
    # write.csv(ABCS_Test, TestOutFile, row.names=FALSE)

} else if (DCP$MODE == "FILE") {
    ABCS_Training <- read.csv(DCP$ABCS_TRAINING_FILE,
                              stringsAsFactors=FALSE)
    ABCS_Test     <- read.csv(DCP$ABCS_TEST_FILE,
                              stringsAsFactors=FALSE)
                              
} else if (DCP$MODE == "CHALLENGE") {
    source(ABCSLOADER)
    ABCS_Training <- loadABCSmaster()
    ABCS_Test     <- read.csv(DCP$ABCS_TEST_FILE,
                              stringsAsFactors=FALSE)
}


# Postcondition: ABCS_Training and ABCS_Test are defined.


# == CREATE FEATURE LISTS ==

# Precondition: ABCS_Training and ABCS_Test are defined.
# EITHER make feature files,
# OR read features from file.

DCP$MAKE_FEATURES <- TRUE

# DCP$DRUG_FEATURE_FILE <- "myDrugFeatures.csv"
# DCP$CELL_FEATURE_FILE <- "myCellFeatures.csv"
# DCP$DRUG_DRUG_FEATURE_FILE <- "myDrugDrugFeatures.csv"
# DCP$DRUG_CELL_FEATURE_FILE <- "myDrugCellFeatures.csv"
# DCP$DRUG_DRUG_CELL_FEATURE_FILE <- "myDrugDrugCellFeatures.csv"

AllFeatures <- list()

if (DCP$MAKE_FEATURES) {

    source(FEATUREMAKER)
    # Featuremaker functions come in five possible flavours:
    # Drug features only, cell features only, drug/drug-,
    # drug/cell-, and drug/drug/cell- combinations.
    # The return values are dataframes:
    #    The first three columns are names; 
    #    Drug A name / Drug B name / Cell line name
    #    These columns are empty if they don't apply but
    #    must be present.
    #    The names are followed by feature-columns.

    AllFeatures$D   <- makeDrugFeatures(ABCS_Training)    
    AllFeatures$C   <- makeCellFeatures(ABCS_Training)
#    AllFeatures$DD  <- makeDrugDrugFeatures(ABCS_Training)     
#    AllFeatures$DC  <- makeDrugCellFeatures(ABCS_Training)     
#    AllFeatures$DDC <- makeDrugDrugCellFeatures(ABCS_Training)     

    # Uncomment if you want to store the results
    # OUT <- sprintf("%s%s_D_Features.csv", DCP$PATH, DCP$PREFIX)
    # write.csv(AllFeatures$D,   OUT, row.names=FALSE)
    #
    # OUT <- sprintf("%s%s_C_Features.csv", DCP$PATH, DCP$PREFIX)
    # write.csv(AllFeatures$C,   OUT, row.names=FALSE)
    #
    # OUT <- sprintf("%s%s_DD_Features.csv", DCP$PATH, DCP$PREFIX)
    # write.csv(AllFeatures$DD,  OUT, row.names=FALSE)
    #
    # OUT <- sprintf("%s%s_DC_Features.csv", DCP$PATH, DCP$PREFIX)
    # write.csv(AllFeatures$DC,  OUT, row.names=FALSE)
    #
    # OUT <- sprintf("%s%s_DDC_Features.csv", DCP$PATH, DCP$PREFIX)
    # write.csv(AllFeatures$DDC, OUT, row.names=FALSE)

} else {
    AllFeatures$D   <- read.csv(DCP$DRUG_FEATURE_FILE, stringsAsFactors=FALSE)
    AllFeatures$C   <- read.csv(DCP$CELL_FEATURE_FILE, stringsAsFactors=FALSE)
#    AllFeatures$DD  <- read.csv(DCP$DRUG_DRUG_FEATURE_FILE, stringsAsFactors=FALSE)
#    AllFeatures$DC  <- read.csv(DCP$DRUG_CELL_FEATURE_FILE, stringsAsFactors=FALSE)
#    AllFeatures$DDC <- read.csv(DCP$DRUG_DRUG_CELL_FEATURE_FILE, stringsAsFactors=FALSE)

}

# Postcondition: Feature dataframes are defined and
#     collected in a list object.


# == COMBINE FEATURE DATA INTO PREDICTION INPUT ==

# Precondition: ABCS_Training set and Feature lists are defined.

# EITHER combine features for each training combination ...
DCP$COMBINE_FEATURES <- TRUE

# OR read combined feature sets from file.
# DCP$COMBI_FEATURE_FILE <- "myCombinedFeatures.csv"

# Note: Obviously, combined features are only relevant for
#       one specific training set.


if (DCP$COMBINE_FEATURES) {

    source(FEATURECOMBINER)
    # Feature combiner functions take as input the ABCS training
    # data and the All_features list.
    # They return a dataframe:
    #    The first three columns are names; 
    #    Drug A name / Drug B name / Cell line name
    #    Column four is a synergy score.
    #    Columns five to seven are synergy categories.
    #    The remaining columns are features.

    CombiFeatures <- makeCombiFeatures(ABCS_Training,
                                       AllFeatures)    

    # Uncomment if you want to store the results
    # OUT <- sprintf("%s%s_CombiFeatures.csv", DCP$PATH, DCP$PREFIX)
    # write.csv(CombiFeatures, OUT, row.names=FALSE)

} else {
    CombiFeatures <- read.csv(DCP$COMBI_FEATURE_FILE, stringsAsFactors=FALSE)
}

# Postcondition: For each training set combination, a
#                feature set has been defined.


# == II: PREDICTIONS ======================================

source(PREDICTOR)

# == TRAIN AND PREDICT TEST DATA ==

PredictionResults <- runPrediction(CombiFeatures,
                                   ABCS_Test)

# runPrediction() takes as input the combined feature set
# and the test-data list.
# It returns a dataframe:
#    The first three columns are names; 
#    Drug A name / Drug B name / Cell line name
#    Column four is the predicted synergy score.
#    Columns five to seven are booleans for synergy
#       category membership.
#    Column eight is the confidence score for the prediction.


# == III: SUBMISSION ======================================

DCP$CHALLENGE_TYPE <- "1A"
# DCP$CHALLENGE_TYPE <- "1B"
# DCP$CHALLENGE_TYPE <- "2"

DCP$PREDICTION_FILE <- sprintf("%s%s_prediction.csv", DCP$PATH, DCP$PREFIX)
DCP$CONFIDENCE_FILE <- sprintf("%s%s_confidence.csv", DCP$PATH, DCP$PREFIX)

source(FORMATTER)

predictionSubmission <- format_prediction(DCP$CHALLENGE_TYPE, PredictionResults)
confidenceSubmission <- format_confidence(DCP$CHALLENGE_TYPE, PredictionResults)

write.csv(predictionSubmission, DCP$PREDICTION_FILE, row.names=FALSE)
write.csv(confidenceSubmission, DCP$CONFIDENCE_FILE, row.names=FALSE)


# == IV: VALIDATION =======================================


# == CALCULATE SCORES ==

source(SCORER)

prediction <- read.csv(DCP$PREDICTION_FILE, stringsAsFactors=FALSE)
confidence <- read.csv(DCP$CONFIDENCE_FILE, stringsAsFactors=FALSE)

scores <- score(DCP$CHALLENGE_TYPE, prediction, confidence)


# == CALCULATE ANALYTICS ==
# TBD

print(scores)



#[END]
