# loadABCSFromMaster.0.1.R
#
# Purpose: Provides functions to load A-B-Cell-Synergy lists
#              from a master source.
#          This file can be loaded as an ABCSLOADER asset in
#              the DREAM_Combi_main.R workflow and is
#              therefore expected to provide a
#              loadABCSmaster() function that returns an
#              ABCS_Master dataframe.
#
#          This file should contain only assets of functions
#             and constants. Sourcing this file must not have
#             any side-effects. Functions should not have side-
#             effects either.
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

loadABCSmaster <- function(IN = "../Challenge Data/Drug Synergy Data/ch1_train_combination_and_monoTherapy.csv") {
    # Return a dataframe with four columns:
    #    Drug A name
    #    Drug B name
    #    Cell line name
    #    Synergy score
    
    master <- read.csv(IN, stringsAsFactors=FALSE)
    master <- master[master[,"QA"] == 1, ] # Exclude all combinations with poor QA scores

    ABCS <- data.frame("A" = master$COMPOUND_A,
                       "B" = master$COMPOUND_B,
                       "C" = master$CELL_LINE,
                       "S" = as.numeric(master$SYNERGY_SCORE),
                       stringsAsFactors = FALSE)
    return(ABCS)
}


# [END]