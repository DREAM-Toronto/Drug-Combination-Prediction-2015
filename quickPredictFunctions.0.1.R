# quickPredictFunctions.0.1.R
#
# Purpose: Provides functions to use our quick predict
#              extrapolation method for synergy score
#              prediction.
#
#          This file can be loaded as a PREDICTOR asset in
#              the DREAM_Combi_main.R workflow and is
#              therefore expected to provide the
#              runPrediction() function that returns 
#              synergy predictions and confidence scores for
#              a test set.
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

runPrediction <- function(FL,       # List of drug and cell feature matrices
                          CF,       # Combined feature set
                          topN = 5, # Use top N best similarities to calculate
                                    # weighted average for similarity score 
                          catProbs = c(0.2, 0.8) # quantile probabilities
                                                 # to assign categories
                          ) {

    # This version of runPrediction() takes as input a list of
    # drug and cell feature matrices and a Combined Feature set.
    
    # All combinations for which synergy scores are reported
    # are used for training. All combinations without
    # scores are to be predicted.
    
    # The function returns a dataframe of predictions only:
    #    The first three columns are names; 
    #    Drug A name / Drug B name / Cell line name
    #    Column four is the predicted synergy score.
    #    Columns five to seven are 0 and 1 values for synergy
    #       category membership.
    #    Column eight is the confidence score for the prediction [0,1].
    
    # FL is expected to contain two objects: 
    #    FL$D is a matrix of drug features,
    #       where the rownames are unique drug names;
    #    FL$C is a matrix of cell features,
    #       where the rownames are unique cell-line names.
    
    # Categories are defined as 1 or 0, depending on whether
    #     the observed synergy is smaller than the lower
    #     of the two quantiles passed in catProbs (col 5),
    #     larger than the higher quantile (col 7), or in- 
    #     between (col 6).
    

    # == extract training set from CF ==
    trainSet <- CF[!(is.na(CF$S)), 1:4]


    # == extract test set from CF ==
    testSet <- CF[is.na(CF$S), 1:7]
    testSet <- data.frame(testSet,
                          "conf" = rep(0, nrow(testSet)),
                          stringsAsFactors = FALSE)

    testSet <- PR$removeUnknown(testSet, FL)
    
    DSM <- PR$makeFeatureSimilarityMatrix(FL$D)
    CSM <- PR$makeFeatureSimilarityMatrix(FL$C)

    testSet <- PR$predictSynergy(trainSet, testSet, DSM, CSM, topN, catProbs)

    return(testSet)
}

# == INTERNAL FUNCTIONS ====================================

PR <- list()

PR$removeUnknown <- function(ABC, FL) {
    # During xValidation testing, it could happen that
    # the test set contains drugs or cells that were
    # not present in the training set. This should not
    # happen in the real submissions, and it should be
    # infrequent enough that we can simply remove the
    # offending row(s) from the test data.

    # If ABC contains drugs or
    # cells that are not in the feature lists, we 
    # collect the row indices.
            
    x <- which(!(ABC$A %in% rownames(FL$D)))
    x <- c(x, which(!(ABC$B %in% rownames(FL$D))))
    x <- c(x, which(!(ABC$C %in% rownames(FL$C))))
    x <- unique(x)
    
    if (length(x) > 0) {
        ABC <- ABC[-(x), ]
    }
    return(ABC)
}


PR$makeFeatureSimilarityMatrix <- function(M) {
    
    N <- nrow(M)
    SM <- matrix(numeric(N * N), ncol=N, nrow=N)
    rownames(SM) <- rownames(M)
    colnames(SM) <- rownames(M)
    for (i in 1:N) {
        for (j in i:N) {
                if (i == j) {
                        SM[i, j] <- 1.0
                } else {
                    cXY <- cor(M[i, ], M[j, ])
                    pXY <- PR$cor2p(cXY)
                SM[i,j] <- SM[j,i] <- pXY
            }
        }
    }
    return(SM)
}


PR$cor2p <- function(C) {
    # Convert a feature correlation coefficient
    # to an estimate of probability that a
    # similar synergy score will be achieved. This
    # particular function is a generalized logistic
    # and entirely empirical with no deeper 
    # theoretical justification.
    
    B <- -11.9
    M <- 0.38
    nu <- 7.26
    
    # x <- seq(-1, 1, by=0.01)
    # plot(x, 1/(1+(exp(B*(x-M))))^nu, type="l", col="#CC2200")
    # abline(v=c(0, 0.5))
    
    return(1/(1+(exp(B*(C-M))))^nu)
}


PR$predictSynergy <- function(train, test, DSM, CSM, topN, catProbs) {
    
    nTrain <- nrow(train)
    nTest <- nrow(test)
    
    simSyn <- data.frame("pSim" = numeric(nTrain),
                         "S"    = numeric(nTrain),
                         stringsAsFactors = FALSE)

    for (iTest in 1:nTest) {
        A_test <- test$A[iTest]
        B_test <- test$B[iTest]
        C_test <- test$C[iTest]
        for (iTrain in 1:nTrain) {
            pA <- DSM[A_test, train$A[iTrain]]  
            pB <- DSM[B_test, train$B[iTrain]]  
            pC <- CSM[C_test, train$C[iTrain]]
            simSyn$pSim[iTrain] <- pA * pB * pC
            simSyn$S[iTrain] <- train$S[iTrain]
        }
        # order descending by pSim
        simSyn <- simSyn[order(simSyn$pSim, decreasing=TRUE), ]
        
        # normalize topN probabilities
        sumSyn <- sum(simSyn$pSim[1:topN])
        normSyn <- simSyn$pSim[1:topN] / sumSyn
        
        # calculate weighted average for topN rows
        fracSyn <- normSyn * simSyn$S[1:topN]
        predSyn <- sum(fracSyn)
        
        # store prediction
        test$S[iTest] <- predSyn

        # approximate confidence as un-normalized
        # sum of probabilites
        test$conf[iTest] <- sumSyn
        
        # convert synergy values to categories for testSet
        # according to quantiles taken from training set
        limits <- quantile(train$S, catProbs, na.rm=TRUE)
        for (i in 1:nrow(test)) {
            if (! is.na(test$S[i])) {
                if (test$S[i] <= limits[1]) {
                        test[i, 5] <- 1
                } else if (test$S[i] >= limits[2]) {
                        test[i, 7] <- 1
                } else {
                        test[i, 6] <- 1                
                }
            }
        }
    }
    return(test)
}


# [END]