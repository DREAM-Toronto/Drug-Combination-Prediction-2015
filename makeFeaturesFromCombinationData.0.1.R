# makeFeaturesFromCombinationData.0.1.R
#
# Purpose: Provides functions to compile unary, binary
#              and ternary features and combinations for
#              all elements provided in an ABCS object.
#
#          This code extracts features from the 
#
#          This file can be loaded as a FEATUREMAKER asset in
#              the DREAM_Combi_main.R workflow and is
#              therefore expected to provide some or all of 
#                 makeDrugFeatures(ABCS)    
#                 makeCellFeatures(ABCS)
#                 makeDrugDrugFeatures(ABCS)     
#                 makeDrugCellFeatures(ABCS)     
#                 makeDrugDrugCellFeatures(ABCS)
#              where ABCS is a dataframe of drug-, drug-,
#              cell-names and synergy scores (for training
#              data), or NA values for test data. 
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


makeDrugFeatures <- function(ABCS,
                             featureSourceFile = "../Challenge Data/Drug Synergy Data/ch1_train_combination_and_monoTherapy.csv",
                             makeUniqueMethod = "independentMeans",
                             imputeMethod = "randomSample",           # set to "NONE" to skip
                             reduceDimMethod = "PCA",                 # set to "NONE" to skip
                             reduceDimTarget = 0.8,                   # cumulative proportion
                             rescaleMethod = "NONE"                   # can be "unitScale"
                             ) {
    
    source <- read.csv(featureSourceFile, stringsAsFactors=FALSE)
    source <- source[source[,"QA"] == 1, ] # Exclude all combinations with poor QA scores

    targetDrugs <- sort(unique(c(ABCS$A, ABCS$B)))
    drugMono <- MF$extractMono(source, "drugs", targetDrugs)
    drugMono <- MF$normalizeConc(drugMono)
    drugMono <- MF$makeUnique(drugMono,
                              method = makeUniqueMethod)
    drugFeatures <- MF$makeDrugFeatureMatrix(drugMono)
    drugFeatures <- MF$imputeMissingTrials(drugFeatures,
                                           method = imputeMethod)
    drugFeatures <- MF$reduceDimensions(drugFeatures,
                                        method = reduceDimMethod,
                                        target = reduceDimTarget)
    drugFeatures <- MF$rescaleFeatures(drugFeatures,
                                       method = rescaleMethod)

    return(drugFeatures)
}


makeCellFeatures <- function(ABCS,
                             featureSourceFile = "../Challenge Data/Drug Synergy Data/ch1_train_combination_and_monoTherapy.csv",
                             makeUniqueMethod = "independentMeans",
                             imputeMethod = "randomSample",           # set to "NONE" to skip
                             reduceDimMethod = "PCA",                 # set to "NONE" to skip
                             reduceDimTarget = 0.8,                   # cumulative proportion
                             rescaleMethod = "NONE"                   # can be "unitScale"
                             ) {
    
    source <- read.csv(featureSourceFile, stringsAsFactors=FALSE)
    source <- source[source[,"QA"] == 1, ] # Exclude all combinations with poor QA scores

    targetCells <- sort(unique(ABCS$C))
    cellMono <- MF$extractMono(source, "cells", targetCells)
    cellMono <- MF$normalizeConc(cellMono)
    cellMono <- MF$makeUnique(cellMono,
                              method = makeUniqueMethod)
    cellFeatures <- MF$makeCellFeatureMatrix(cellMono)
    cellFeatures <- MF$imputeMissingTrials(cellFeatures,
                                           method = imputeMethod)
    cellFeatures <- MF$reduceDimensions(cellFeatures,
                                        method = reduceDimMethod,
                                        target = reduceDimTarget)
    cellFeatures <- MF$rescaleFeatures(cellFeatures,
                                       method = rescaleMethod)

    return(cellFeatures)
}



# == INTERNAL FUNCTIONS ====================================

MF <- list()

MF$extractMono <- function(source, type, targets) {
	
	N <- 2 * nrow(source)
	mono <- data.frame("D" = character(N),
	                   "C" = character(N),
	                   "max" = numeric(N),
	                   "IC50" = numeric(N),
	                   "H" = numeric(N),
	                   "Einf" = numeric(N),
	                   stringsAsFactors = FALSE)
	j <- 0
	for (i in 1:nrow(source)) {
		if (type == "drugs" & (source$COMPOUND_A[i] %in% targets)) {
			j <- j+1
			mono$D[j]    <- source$COMPOUND_A[i]
			mono$C[j]    <- source$CELL_LINE[i]
			mono$max[j]  <- source$MAX_CONC_A[i]
			mono$IC50[j] <- source$IC50_A[i]
			mono$H[j]    <- source$H_A[i]
			mono$Einf[j] <- source$Einf_A[i]
		}
		if (type == "drugs" & (source$COMPOUND_B[i] %in% targets)) {
			j <- j+1
			mono$D[j]    <- source$COMPOUND_B[i]
			mono$C[j]    <- source$CELL_LINE[i]
			mono$max[j]  <- source$MAX_CONC_B[i]
			mono$IC50[j] <- source$IC50_B[i]
			mono$H[j]    <- source$H_B[i]
			mono$Einf[j] <- source$Einf_B[i]
		}
		if (type == "cells" & (source$CELL_LINE[i] %in% targets)) {
			j <- j+1
			mono$D[j]    <- source$COMPOUND_A[i]
			mono$C[j]    <- source$CELL_LINE[i]
			mono$max[j]  <- source$MAX_CONC_A[i]
			mono$IC50[j] <- source$IC50_A[i]
			mono$H[j]    <- source$H_A[i]
			mono$Einf[j] <- source$Einf_A[i]
			j <- j+1
			mono$D[j]    <- source$COMPOUND_B[i]
			mono$C[j]    <- source$CELL_LINE[i]
			mono$max[j]  <- source$MAX_CONC_B[i]
			mono$IC50[j] <- source$IC50_B[i]
			mono$H[j]    <- source$H_B[i]
			mono$Einf[j] <- source$Einf_B[i]
		}
	}
	return(mono[1:j, ])
}


MF$normalizeConc <- function(data) {
	# normalize drug concentrations from 0 to maxConcentration = 1
	
	drugs <- sort(unique(data$D))
	maxC <- rep(0, length(drugs))
	names(maxC) <- drugs
	
	for (i in 1:nrow(data)) {
		maxC[data$D[i]] <- max(maxC[data$D[i]], data$max[i])
	}

	for (i in 1:nrow(data)) {
		data$IC50[i] <- data$IC50[i] / maxC[data$D[i]]
		data$H[i]    <- data$H[i]    * maxC[data$D[i]]
	}
	
	return(data)
}


MF$makeUnique <- function(data, method) {
	# Return a data set with unique DC value combinations.
	# Handle multiple observations according to "method".
	
	combi <- unique(data.frame("D"    = data$D,
	                           "C"    = data$C,
	                           "IC50" = numeric(nrow(data)),
	                           "H"    = numeric(nrow(data)),
	                           "Einf" = numeric(nrow(data)),
	                           stringsAsFactors = FALSE))
	
	for (i in 1:nrow(combi)) {
		rows <- data[data$D == combi$D[i] & data$C == combi$C[i], ]
		if (method == "independentMeans") {
			combi$IC50[i] <- mean(rows$IC50)
			combi$H[i]    <- mean(rows$H)
			combi$Einf[i] <- mean(rows$Einf)
		}
	}
	rownames(combi) <- 1:nrow(combi)
	
	return(combi)
}


MF$makeDrugFeatureMatrix <- function(drugMono) {
	
	drugs <- unique(drugMono$D)
	cells <- unique(drugMono$C)
	DFM <- matrix(rep(NA, length(drugs) * length(cells) * 3),
	              nrow=length(drugs))
	rownames(DFM) <- drugs
	
	for (i in 1:length(drugs)) {
		for (j in 1:length(cells)) {
			ii <- which(drugMono$D == drugs[i] &
			            drugMono$C == cells[j])
			if(length(ii) > 0) {  # combination has been observed ...
    			jj <- ((j*3)-2)
	     		DFM[i, jj]   <- drugMono$IC50[ii]
		     	DFM[i, jj+1] <- drugMono$H[ii]
			    DFM[i, jj+2] <- drugMono$Einf[ii]				
			}
		}
	}
	return(DFM)
}


MF$makeCellFeatureMatrix <- function(cellMono) {
	
	drugs <- unique(cellMono$D)
	cells <- unique(cellMono$C)
	CFM <- matrix(rep(NA, length(drugs) * length(cells) * 3),
	              nrow=length(cells))
	rownames(CFM) <- cells
	
	for (i in 1:length(cells)) {
		for (j in 1:length(drugs)) {
			ii <- which(cellMono$C == cells[i] &
			            cellMono$D == drugs[j])
			if(length(ii) > 0) {  # combination has been observed ...
    			jj <- ((j*3)-2)
	     		CFM[i, jj]   <- cellMono$IC50[ii]
		     	CFM[i, jj+1] <- cellMono$H[ii]
			    CFM[i, jj+2] <- cellMono$Einf[ii]				
			}
		}
	}
	return(CFM)
}


MF$imputeMissingTrials <- function(FM, method) {

    if (method == "randomSample") {
		base <- seq(1, ncol(FM), by = 3)
		for (i in 1:nrow(FM)) {
			rowIC50 <- FM[i, base]
			rowIC50 <- rowIC50[!is.na(rowIC50)]
		
			rowH <- FM[i, base + 1]
			rowH <- rowH[!is.na(rowH)]
		
			rowEinf <- FM[i, base + 2]
			rowEinf <- rowEinf[!is.na(rowEinf)]
		
		    for (j in base) {
			    	if (is.na(FM[i, j])) {
			    		FM[i, j]     <- sample(rowIC50, 1)
			    		FM[i, j + 1] <- sample(rowH,    1)
			    		FM[i, j + 2] <- sample(rowEinf, 1)
			    	}
		    }
		}
        return(FM)
        
    } else if (method == "NONE") {
    	return(FM)
    	
    } else {
    	stop(paste("PANIC: unsupported method", method, "requested."))
    }
}


MF$reduceDimensions <- function(M, method, target) {

    if (method == "PCA") {
		for (i in 1:ncol(M)) {
			M[ , i] <- scale(M[ , i])
		}
		pcaM <- prcomp(M)
		# plot(pcaM)
		# summary(pcaM)
        cols <- summary(pcaM)$importance[3, ] < target
        return(pcaM$x[ , cols])
        
    } else if (method == "NONE") {
    	return(M)
    	
    } else {
    	stop(paste("PANIC: unsupported method", method, "requested."))
    }
}


MF$rescaleFeatures <- function(M, method) {

    if (method == "unitScale") {
    	# Note: unit scale may be required for ML input,
    	# but may be misleading for the general case
    	# since it rescales the PCA-derived loadings. 
    	
    	for (i in 1:ncol(M)) {
	        x <- M[ , i]                   # fetch a column
	        x <- (x / (max(x) - min(x)))   # re-scale
	        x <- x - min(x)                # re-center
	        M[ , i] <- x                   # re-assign
	    }
        return(M)
        
    } else if (method == "NONE") {
    	return(M)
    	
    } else {
    	stop(paste("PANIC: unsupported method", method, "requested."))
    }
}



# [END]