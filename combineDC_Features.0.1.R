# combineDC_Features.0.1.R
#
# Purpose: Provides functions to combine Drug and Cell features
#              as input for machine learning.
#
#          This file can be loaded as a FEATURECOMBINER asset in
#              the DREAM_Combi_main.R workflow and is
#              therefore expected to provide the
#              makeCombiFeatures() function that returns a
#              Feature set.
#
#          This file should contain only assets of functions
#             and constants. Sourcing this file must not have
#             any side-effects. Functions should not have side-
#             effects either.
#
# Version: 0.1           
#
# Date:    Jan 23 2016
# Author:  Boris and DREAM team UofT
#
# V 0.1    First code
#
# TODO:    
#
# ==========================================================

makeCombiFeatures <- function(CF,  # ABCS input dataframe
                              FL,  # List of feature matrices
                              catProbs = c(0.2, 0.8) # quantile probabilities
                                                     # to assign categories
                              ) {

    # Returns a dataframe
    #    The first three columns are names; 
    #    Drug A name / Drug B name / Cell line name
    #    Column four is a synergy score.
    #    Columns five to seven are synergy categories.
    #    The remaining columns are features.
    
    # Unobserved features or combinations are coded as NA.
    
    # This version expects FL to contain two objects: 
    #    FL$D is a matrix of drug features,
    #       where the rownames are unique drug names;
    #    FL$C is a matrix of cell features,
    #       where the rownames are unique cell-line names.
    
    # Categories are defined as 1 or 0, depending on whether
    #     the observed synergy is smaller than the lower
    #     of the two quantiles passed in catProbs (col 5),
    #     larger than the higher quantile (col 7), or in- 
    #     between (col 6).
    
    # == add required number of empty columns==
    CF <- data.frame(CF,
                     matrix(0, nrow=nrow(CF), ncol=3),
                     matrix(as.numeric(NA), nrow=nrow(CF), ncol=ncol(FL$D)), # drugs
                     matrix(as.numeric(NA), nrow=nrow(CF), ncol=ncol(FL$D)), # drugs
                     matrix(as.numeric(NA), nrow=nrow(CF), ncol=ncol(FL$C)), # cells
                     stringsAsFactors = FALSE)

    # == process CF and add all features ==
    
    first <- 8
    last <- ncol(CF)
    
    for (i in 1:nrow(CF)) {
    	    CF[i, first:last] <- c(FL$D[CF$A[i], ], FL$D[CF$B[i], ], FL$C[CF$C[i], ])
    }

    # == define categories ==
    
    limits <- quantile(CF$S, catProbs, na.rm=TRUE)
    for (i in 1:nrow(CF)) {
    	    if (! is.na(CF$S[i])) {
	    	    if (CF$S[i] <= limits[1]) {
	    	    	    CF[i, 5] <- 1
	    	    } else if (CF$S[i] >= limits[2]) {
	    	    	    CF[i, 7] <- 1
	    	    } else {
	    	    	    CF[i, 6] <- 1    	    	
	    	    }
    	    }
    	}
    return(CF)
}


# [END]