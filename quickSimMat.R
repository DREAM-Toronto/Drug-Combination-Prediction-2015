# quickSimMat.R
# Quick calculation of similarity matrices for drugs and cells
#
# Puts DSM and CSM similarity matrices into the environment
#
# To Do:   This actually will make drugs / cells with fewer
#             observations more similar to each other ...
#             calculate correlations only for pairs that
#             are not both equal to column average (i.e.
#             unobserved).
#
# Version: 0.3           
#
# Date:    Nov 9 2015
# Author:  Boris and DREAM team UofT
#          
# V 0.3    Refactor code for easier alternative handling
#              of correlation calculations.
#          Change correlation calculations to ignore pairs
#              of unobserved values, to prevent pairs with
#              many missing features to appear more similar.
#              This improves the prediction quality by
#              about 25%.
# V 0.2.1  Add progress messages.
# V 0.2    Maintenance and refactoring.
#          Bugfix in calculating CSM
# V 0.1    First code
# ==========================================================

# This file needs to be source()'d from predictSynergy.R 

# == FUNCTIONS =============================================
#

uniqueCor <- function(x, y) {
	# x and y are feature vectors of triplets of 
	# (IC50, H, Einf) for each drug or cell line.
	#   - Remove all triplets in which both x and y
	#     have NAs; 
	#   - replace the NAs with averages;
	#   - calculate correlation;
	#   - convert correlation to probability.
	# Output is the probability of similarity.
    obs <- !(is.na(x) & is.na(y))
    x <- x[obs]
    y <- y[obs]
    x <- na2means(x) 
    y <- na2means(y)
    cXY <- cor(x, y)
    pXY <- cor2p(cXY)
    return(pXY) 
}

na2means <- function(x) {
	# replace all NA triplets in vector with
	# averages
    tri <- seq(1, length(x), by=3)
    xMeans <- c(mean(x[tri],   na.rm=TRUE),
                mean(x[tri+1], na.rm=TRUE),
                mean(x[tri+2], na.rm=TRUE))
	if (is.na(xMeans[1])) {
		stop("PANIC: feature means is NA in vector.")
    }
    z <- rep(xMeans, length(tri))
    sel <- is.na(x)
    x[sel] <- z[sel]
    return(x)
}

cor2p <- function(cor) {
	# crude manual approximation to a CDF that converts
	# a correlation coefficient to a probability (of
	# having similar synergy score)
	if(cor < 0.1) { return(0)}
	if(cor < 0.2) { return(0.01)}
	if(cor < 0.3) { return(0.02)}
	if(cor < 0.4) { return(0.05)}
	if(cor < 0.5) { return(0.1)}
	if(cor < 0.6) { return(0.2)}
	if(cor < 0.7) { return(0.3)}
	if(cor < 0.8) { return(0.7)}
	if(cor < 0.9) { return(0.9)}
                    return(1.0)
}

# == MAIN ==================================================
#

if (VERBOSE) {cat(paste("    preparing features ..."))}

trainingData <- read.csv(TRAINING_SET_FILE,
                         header = TRUE,
                         row.names = NULL,
                         stringsAsFactors = FALSE)
# head(trainingData)

# == PREPARE FEATURE MATRIX =======

# Read all compounds and cells
nTrain <- nrow(trainingData)
cells <- character(nTrain)
drugs <- character(nTrain * 2)
for (i in 1:nTrain) {
	cells[i]        <- trainingData[i, "CELL_LINE"]
	drugs[i]        <- trainingData[i, "COMPOUND_A"]
	drugs[nTrain+i] <- trainingData[i, "COMPOUND_B"]
}

# Make compounds and cells unique
cells <- sort(unique(cells))
drugs <- sort(unique(drugs))

nCells <- length(cells)
nDrugs <- length(drugs)

# Make a matrix of entries for monotherapies - we store the
# triplets for Einf, IC50 and H as adjacent columns.
data <- matrix(NA, nrow=nDrugs, ncol=nCells * 3)

# Set the row names to drugs
rownames(data) <- drugs

# vector of indices of the first columns of each parameter triplet
colIC50 <- (((1:nCells)-1) * 3) + 1

# Set the column names to <cell>.<param>
cn <- character(nCells * 3)
for (i in 1:nCells) {
	ii <- colIC50[i]
	cn[ii] <- cells[i] # this column holds IC_50 values; use only
	                   # the cell_line name, so we don't need to 
	                   # regex the name out of the colname downstream
    cn[ii+1] <- paste(cells[i], ".H", sep="")
    cn[ii+2] <- paste(cells[i], ".Einf", sep="")
}
colnames(data) <- cn
# head(data)


# Iterate over all compound / cell monotherapies.
# If there is more than one monotherapy for a cell
# then average the DRC parameters for all of them.
for (id in 1:nDrugs) {
	for (ic in 1:nCells) {
        # for each monotherapy, collect all three values
		mono_A <- trainingData[trainingData["CELL_LINE"] == cells[ic] &
		                       trainingData["COMPOUND_A"] == drugs[id],
		                       c("IC50_A", "H_A", "Einf_A")]
		# change colnames to prevent names mismatch when rowbind()'ing
		# data for compound B                     
		colnames(mono_A) <- c("IC50_B", "H_B", "Einf_B") 

		mono_B <- trainingData[trainingData["CELL_LINE"] == cells[ic] &
		                       trainingData["COMPOUND_B"] == drugs[id],
		                       c("IC50_B", "H_B", "Einf_B")]
        mono <- rbind(mono_A, mono_B)
	    if (nrow(mono) > 0) {
	        # ... average the monotherapies, (in case there
	        #     was more than one),		
	        IC50 <- mean(mono[,1])
	        H    <- mean(mono[,2])
	        Einf <- mean(mono[,3])
            # ... and write the DRC parameter into the data matrix
            data[drugs[id], cells[ic]] <- IC50
            data[drugs[id], paste(cells[ic], ".H", sep="")] <- H
            data[drugs[id], paste(cells[ic], ".Einf", sep="")] <- Einf  
		} 
	}
}
# head(data)


# == COMPILE DSM DRUG SIMILARITY MATRIX =======

if (VERBOSE) {cat(paste(" compiling DSM ..."))}

# initialize the DSM
DSM <- matrix(numeric(nDrugs * nDrugs), ncol=nDrugs, nrow=nDrugs)
rownames(DSM) <- drugs
colnames(DSM) <- drugs

# calculate correlations for each pair of feature vectors
# and store in DSM
for (i in 1:nDrugs) {
	for (j in i:nDrugs) {
		DSM[i,j] <- DSM[j,i] <- uniqueCor(data[i, ], data[j, ])
	}
}
# head(DSM)	


# == COMPILE CSM CELL SIMILARITY MATRIX =======

if (VERBOSE) {cat(paste(" compiling CSM ..."))}

# initialize the CSM
CSM <- matrix(numeric(nCells * nCells), ncol=nCells, nrow=nCells)
rownames(CSM) <- cells
colnames(CSM) <- cells

# collapse the column triplets into feature vectors
fv <- matrix(numeric(nCells * nDrugs * 3),
             nrow = nCells,
             ncol = (nDrugs * 3))             
for (i in 1:nCells) {
	j <- colIC50[i]
	# each vector is joined from the parameter triplets
	# stored for that cell, for each drug
	fv[i, ] <- as.vector(t(data[ , j:(j+2)]))
}
# head(fv)

# calculate correlations for each pair of feature vectors
# and store in CSM
for (i in 1:nCells) {
	for (j in i:nCells) {
		CSM[i,j] <- CSM[j,i] <- uniqueCor(fv[i, ], fv[j, ])
	}
}
# head(CSM)	

if (VERBOSE) {cat(paste(" Done\n"))}
# Done: DSM and CSM have been computed


# [END]
