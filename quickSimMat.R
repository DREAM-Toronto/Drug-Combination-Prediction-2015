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
# Version: 0.2           
#
# Date:    Nov 7 2015
# Author:  Boris and DREAM team UofT
#          
# V 0.2    Maintenance and refactoring.
#          Bugfix in calculating CSM
# V 0.1    First code
# ==========================================================

# This file needs to be source()'d from predictSynergy.R 

# == FUNCTIONS =============================================
#

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

d2 <- data  # make a copy that we can modify
            # by replacing NAs with
            # a value estimate

# replace all NAs with row-averages
for (i in 1:nrow(d2)) {
	# calculate row averages
	mIC50 <- mean(d2[i, colIC50  ],   na.rm=TRUE)
	mH    <- mean(d2[i, colIC50 + 1], na.rm=TRUE)
	mEinf <- mean(d2[i, colIC50 + 2], na.rm=TRUE)
	
	if (is.nan(mIC50)) {
		stop(sprintf("PANIC: no data at all for compound in row %d ", i))
	}
	
	for (j in colIC50) {
		# replace values if IC50 is NA
		if (is.na(d2[i, j])) {
			d2[i, j]   <- mIC50
			d2[i, j+1] <- mH
			d2[i, j+2] <- mEinf
		}
	}
}
# head(d2)

# compile the DSM
DSM <- matrix(numeric(nDrugs * nDrugs), ncol=nDrugs, nrow=nDrugs)
rownames(DSM) <- drugs
colnames(DSM) <- drugs

for (i in 1:nDrugs) {
	for (j in i:nDrugs) {
		# similarity is cor2p() of the correlation
		# coefficient of rows i and j
		DSM[i,j] <- DSM[j,i] <- cor2p(cor(d2[i, ], d2[j, ]))
	}
}
	
# head(DSM)	



# == COMPILE CSM CELL SIMILARITY MATRIX =======

# replace all NAs in "data" with column-averages
for (i in colIC50) {
	mIC50 <- mean(data[ , i  ], na.rm=TRUE)
	mH    <- mean(data[ , i+1], na.rm=TRUE)
	mEinf <- mean(data[ , i+2], na.rm=TRUE)

	if (is.nan(mIC50)) {
		stop(sprintf("PANIC: no data at all for cell line in column %d ", i))
	}

	for (j in 1:nrow(data)) {
		if (is.na(data[j, i])) {
			data[j, i]   <- mIC50
			data[j, i+1] <- mH
			data[j, i+2] <- mEinf
		}
	}
}
# head(data)

# initialize the CSM
CSM <- matrix(numeric(nCells * nCells), ncol=nCells, nrow=nCells)
rownames(CSM) <- cells
colnames(CSM) <- cells

# Note: below is erroneous code which was used in the
# round 1 submission. The column-triplets are not 
# correctly iterated over. It is kept here 
# until the code has been submitted to git so that we have
# a documentation of the error, and can experiment with
# its consequences.
# for (i in 1:n) {
	# cell_i <- c(data[ , i], data[ , i+1], data[ , i+2])
	# for (j in 1:n) {
	    # cell_j <- c(data[ , j], data[ , j+1], data[ , j+2])
		# CSM[i,j] <- cor2p(cor(cell_i, cell_j))
	# }
# }

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
		CSM[i,j] <- CSM[j,i] <- cor2p(cor(fv[i, ], fv[j, ]))
	}
}

# head(CSM)	

# Done: DSM and CSM have been computed


# [END]
