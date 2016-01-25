# formatSubmissionFunctions.0.1.R
#
# Purpose: Provides functions to load format a prediction
#              dataframe into a valid DREAM challenge submission.
#
#          This file can be loaded as a FORMATTER asset in
#              the DREAM_Combi_main.R workflow and is
#              therefore expected to provide a
#              format_prediction() function and a
#              format_confidence() function.
#
#          This file should contain only assets of functions
#             and constants. Sourcing this file must not have
#             any side-effects. Functions should not have side-
#             effects either.
#
# Version: 0.1           
#
# Date:    Jan 24 2016
# Author:  Boris and DREAM team UofT
#
# V 0.1    First code
#
# TODO:    
#
# ==========================================================

format_prediction <- function(type, pred) {
    
    # From https://www.synapse.org/#!Synapse:syn4231880/wiki/235659
    # Constraints:
    #    CELL_LINE column contains cell line identifier. Those IDs
    #    are the normalized cell line name, which is controlled vocabulary.
    #    COMBINATION_ID column contains the combination identifier, 
    #    which consist out of both drug names separated by a dot. 
    #    Please note that drug names are controlled vocabulary too 
    #    and alphabetically sorted.
    #    PREDICTION column contains scalar values of the predicted
    #    synergy, where positive values indicate synergy, values around
    #    zero additivity and negative values infer antagonism.
    #    NA's (i.e. null) predictions are not accepted.
    
    
    if (type == "1A" | type == "1B") {

        N <- nrow(pred)
        sub <- data.frame("CELL_LINE" = character(N),
                          "COMBINATION_ID"  = character(N),
                          "PREDICTION" = numeric(N),
                          stringsAsFactors = FALSE)
        for (i in 1:N) {
                sub$CELL_LINE[i]      <- pred$C[i]
                sub$COMBINATION_ID[i] <- sprintf("%s.%s", pred$A[i], pred$B[i])
                sub$PREDICTION[i]     <- pred$S[i]                
        }
        
    } else {
        stop("PANIC: type \"", type, "\" not yet supported.", sep="")
    }
    
    sub <- sub[order(sub$COMBINATION_ID), ]
    return(sub)

}

format_confidence <- function(type,              # Challenge type
                              pred,              # Prediction input
                              method="quantile", # Use method to calculate
                                                 # single confidence 
                                                 # from multiple combination IDs
                              pQuant=0.8         # Return value of pQuant quantile
                              ) {

    # From https://www.synapse.org/#!Synapse:syn4231880/wiki/235659
    # Constraints:
    #   - Confidence is in a range from 0 to 1, where high values
    #     correspond to high confidence in that particular combination.
    #   - Confidence value cannot be null or NA.
    #   - Only valid COMBINATION_ID's, meaning controlled vocabulary
    #     drug names separated by dot and alphabetically ordered.
    #   - List must contain all combinations from challenge 1.
    #   - No duplicated COMBINATION_ID's allowed!

    if (type == "1A" | type == "1B") {

        N <- nrow(pred)
        conf <- data.frame("A" = character(N),
                           "B" = character(N),
                           "COMBINATION_ID"  = character(N),
                           "CONFIDENCE" = numeric(N),
                           stringsAsFactors = FALSE)

        for (i in 1:N) {
                conf$A[i]              <- pred$A[i]
                conf$B[i]              <- pred$B[i]
                conf$COMBINATION_ID[i] <- sprintf("%s.%s", pred$A[i], pred$B[i])
        }
        conf <- unique(conf)
        conf <- conf[order(conf$COMBINATION_ID), ]
        
        for (i in 1:nrow(conf)) {
            pVals <- pred$conf[pred$A == conf$A[i] & pred$B == conf$B[i]]
            if (method == "quantile") {
                    conf$CONFIDENCE[i] <- quantile(pVals, probs = pQuant)               
            } else {
                    stop("PANIC: method \"", method, "\" not supported.", sep="")
            }
        }
        # rank and normalize
        conf$CONFIDENCE <- rank(conf$CONFIDENCE)/nrow(conf)

        # drop A and B columns
        conf <- conf[ , c("COMBINATION_ID", "CONFIDENCE")]
        
    } else {
        stop("PANIC: type \"", type, "\" not yet supported.", sep="")
    }
    return(conf)
}


# [END]