# DREAM_Combi_main.R
#
# Purpose: Configure and drive a prediction run for the Drug
#             Combination Prediction Challenge. This is a 
#             framework to tie together our code assets.
#
# Version: 0.1           
#
# Date:    Jan 17 2016
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
DCP$PREFIX <- "01.test"   # filename prefix for this run

# == HELPERS ==
source(DREAM_Combi_HelperFunctions.R) 

# == INITIALIZE ============================================


# == I: FEATURES ===========================================


# == DEFINE TRAINING DATA AND TEST DATA ==

# No precondition.
# EITHER make a random training and test set
# OR read ABCS training and test sets from file.

DCP$MAKE_ABCS_TT <- TRUE
DCP$HOLDOUT_FRACTION <- 0.33


# DCP$ABCS_TRAINING_FILE <- "myAwesomeTrainingSet.csv"
# DCP$ABCS_TEST_FILE     <- "myAwesomeTestSet.csv"

if (DCP$MAKE_ABCS_TT) {
    ABCSLOADER <- "loadABCSFromMaster.0.1.R"

    source(ABCSLOADER)

    # Load A-B-C-Scores
    # The return value of loadABCSmaster() is a dataframe
    # with four columns: 
    #    Drug A name
    #    Drug B name
    #    Cell line name
    #    Synergy score

    # set.seed(112358)   # uncomment for reproducible runs
                         # for debugging purposes only!
    ABCSmaster <- loadABCSmaster(...)
	nHold <- round(nrow(ABCSmaster) * DCP$HOLDOUT_FRACTION) # number of holdouts
	iHold <- sample(1:nrow(ABCSmaster), nHold)   # random row index of holdouts
	
	ABCS_Training <- master[-(iHold), ]
	ABCS_Test  <- master[iHold, ]

    # Uncomment if you want to store the results
	# write.csv(ABCS_Training, "tmp.csv", row.names=FALSE)
	# write.csv(ABCS_Test, "tmp.csv", row.names=FALSE)

} else {
	ABCS_Training <- read.csv(DCP$ABCS_TRAINING_FILE,
	                          stringsAsFactors=FALSE)
	ABCS_Test     <- read.csv(DCP$ABCS_TEST_FILE,
	                          stringsAsFactors=FALSE)
}


# Postcondition: ABCS_Training and ABCS_Test are defined.


# == CREATE FEATURE LISTS ==

# Precondition: ABCS_Training and ABCS_Test are defined.
# EITHER make feature files,
# OR read features from file.

DCP$MAKE_FEATURES <- TRUE
DCP$MAKE_DRUG_FEATURES <- TRUE
DCP$MAKE_CELL_FEATURES <- TRUE
DCP$MAKE_DRUG_DRUG_FEATURES <- TRUE
DCP$MAKE_DRUG_CELL_FEATURES <- TRUE
DCP$MAKE_DRUG_DRUG_CELL_FEATURES <- TRUE

# DCP$DRUG_FEATURE_FILE <- "myDrugFeatures.csv"
# DCP$CELL_FEATURE_FILE <- "myCellFeatures.csv"
# DCP$DRUG_DRUG_FEATURE_FILE <- "myDrugDrugFeatures.csv"
# DCP$DRUG_CELL_FEATURE_FILE <- "myDrugCellFeatures.csv"
# DCP$DRUG_DRUG_CELL_FEATURE_FILE <- "myDrugDrugCellFeatures.csv"

AllFeatures <- list()

if (DCP$MAKE_FEATURES) {
    FEATUREMAKER <- "makeFeaturesFromMonotherapy.0.1.R"
    source(FEATUREMAKER)
    # Featuremaker functions come in five possible flavours:
    # Drug features only, cell features only, drug/drug-,
    # drug/cell-, and drug/drug/cell- combinations.
    # The return values are dataframes:
    #    The first three columns are names; 
    #    Drug A name / Drug B name / Cell line name
    #    These columns are empty if they don't apply but
    #    must be present.
    #    The names are followed by feature-columns

    	AllFeatures$D   <- makeDrugFeatures(ABCS_Training, ...)    
    	AllFeatures$C   <- makeCellFeatures(ABCS_Training, ...)
#    AllFeatures$DD  <- makeDrugDrugFeatures(ABCS_Training, ...)     
#    AllFeatures$DC  <- makeDrugCellFeatures(ABCS_Training, ...)     
#    AllFeatures$DDC <- makeDrugDrugCellFeatures(ABCS_Training, ...)     

    # Uncomment if you want to store the results
	# write.csv(AllFeatures$D,   "tmp.csv", row.names=FALSE)
	# write.csv(AllFeatures$C,   "tmp.csv", row.names=FALSE)
	# write.csv(AllFeatures$DD,  "tmp.csv", row.names=FALSE)
	# write.csv(AllFeatures$DC,  "tmp.csv", row.names=FALSE)
	# write.csv(AllFeatures$DDC, "tmp.csv", row.names=FALSE)

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
# EITHER combine features for each training combination,
# OR read combined feature sets from file.
# Note: Obviously, combined features are only relevant for
#       one specific training set.

DCP$COMBINE_FEATURES <- TRUE
# DCP$COMBI_FEATURE_FILE <- "myCombinedFeatures.csv"

if (DCP$COMBINE_FEATURES) {
    FEATURECOMBINER <- "combineDC_Features.0.1.R"
    source(FEATURECOMBINER)
    # Feature combiner functions take as input the ABCS training
    # data and the All_features list.
    # They return a dataframe:
    #    The first three columns are names; 
    #    Drug A name / Drug B name / Cell line name
    #    Column four is a synergy score.
    #    The remaining columns are features.

    	CombiFeatures <- makeCombiFeatures(ABCS_Training,
    	                              AllFeatures,
    	                              ...)    

    # Uncomment if you want to store the results
	# write.csv(CombiFeatures, "tmp.csv", row.names=FALSE)

} else {
    	CombiFeatures   <- read.csv(DCP$COMBI_FEATURE_FILE, stringsAsFactors=FALSE)
}

# Postcondition: For each training set combination, a
#               feature set has been defined.


# == II: PREDICTIONS ======================================



# == TRAIN ALGORITHM ==



# == PREDICT TEST DATA ==


# == III: SUBMISSION ======================================


# == PREPARE PREDICTION SUBMISSION ==


# == PREPARE CONFIDENCE SCORE SUBMISSION ==


# == IV: VALIDATION =======================================


# == CALCULATE SCORES ==


# == CALCULATE ANALYTICS ==




#[END]
