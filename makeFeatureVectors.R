# makeFeatureVectors.R
# Produces feature vector and synergy data files from
# training data.
#
# Version: 0.2           
#
# Date:    Nov 12 2015
# Author:  Boris and DREAM Toronto
#
# V 0.2    Normalize drug concentrations;
#          Remove rows where QA != 1;
#          Add second output file in which NAs are
#              converted to row-averages.
# V 0.1    First code
# ==========================================================

# == CONSTANTS =============================================
#

INPUT_FILE <- "../Challenge Data/Drug Synergy Data/ch1_train_combination_and_monoTherapy.csv"
OUT_NA_FILE <- "drugFeaturesNA.csv"  # contains NA for unobserved
OUT_AV_FILE <- "drugFeaturesAv.csv"  # contains row averages for unobserved
SYNERGY_FILE <- "drugSynergies.csv"

# == FUNCTIONS =============================================
#

# == MAIN ==================================================
#

# == READ AND PREPARE DATA
trainingData <- read.csv(INPUT_FILE, stringsAsFactors=FALSE)
trainingData <- trainingData[trainingData[ , "QA"] == 1, ]


# Extract names of compounds and cells
nTrain <- nrow(trainingData)
cells <-   trainingData[ , "CELL_LINE"]
drugs <- c(trainingData[ , "COMPOUND_A"],
           trainingData[ , "COMPOUND_B"])


# Make list of unique drugs and cells
cells <- sort(unique(cells))
drugs <- sort(unique(drugs))

nCells <- length(cells)
nDrugs <- length(drugs)


# Normalize concentrations: divide IC50 values
# by MAX_CONC
trainingData[ , "IC50_A"] <- trainingData[ , "IC50_A"] / trainingData[ , "MAX_CONC_A"]
trainingData[ , "IC50_B"] <- trainingData[ , "IC50_B"] / trainingData[ , "MAX_CONC_B"]
# head(trainingData)


# Make a matrix of entries for monotherapies - we store the
# triplets for Einf, IC50 and H as adjacent columns.
data <- matrix(NA, nrow=nDrugs, ncol=nCells * 3)

# Set the row names to drugs
rownames(data) <- drugs

# vector of indices of the first columns of each parameter triplet
colIC50 <- seq(1, nCells*3, by = 3)

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


# == EXTRACT FEATURE VECTORS
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
	        IC50 <- mean(mono[ , "IC50_B"])
	        H    <- mean(mono[ , "H_B"])
	        Einf <- mean(mono[ , "Einf_B"])
	        
            # ... and write the DRC parameter into the data matrix
            data[drugs[id], cells[ic]] <- IC50
            data[drugs[id], paste(cells[ic], ".H", sep="")] <- H
            data[drugs[id], paste(cells[ic], ".Einf", sep="")] <- Einf  
		} 
	}
}
# head(data)


# == EXTRACT SYNERGY DATA
syn <- aggregate(SYNERGY_SCORE ~ COMBINATION_ID, trainingData, median)
syn <- data.frame(COMBINATION_ID = syn$COMBINATION_ID, 
                  COMPOUND_A = "",
                  COMPOUND_B = "",
                  SYNERGY_SCORE = syn$SYNERGY_SCORE,
                  stringsAsFactors = FALSE)
# head(syn)

for (i in 1:nrow(syn)) {
	# collect compound A and B names
	row <- (trainingData[trainingData$COMBINATION_ID == syn$COMBINATION_ID[i], ])[1, ]
	syn[i, "COMPOUND_A"] <- row["COMPOUND_A"] 
	syn[i, "COMPOUND_B"] <- row["COMPOUND_B"] 
}
# head(syn)

# Check whether all rows in the feature table is a drug
# that has been observed in a combination at least once.
observedCombDrugs <- unique(c(syn$COMPOUND_A, syn$COMPOUND_B))
v <- c()
for (i in 1:nrow(data)) {
	v[i] <- !any(observedCombDrugs %in% rownames(data)[i])
}
if (sum(v) != 0) {
	stop("unobserved compound in data")
}


# Make a copy of data to replace NA values 
# with row averages

data2 <- data
for (i in 1:nrow(data2)) {
	rowMeans <- c(mean(data2[i, colIC50],   na.rm=TRUE),
                  mean(data2[i, colIC50+1], na.rm=TRUE),
                  mean(data2[i, colIC50+2], na.rm=TRUE))
    for (j in colIC50) {
    	if (is.na(data2[i, j])) {
    		data2[i, j:(j+2)] <- rowMeans
    	}
    }
}
# head(data2)

# == WRITE FILES
write.csv(data,  OUT_NA_FILE, row.names=TRUE)
write.csv(data2, OUT_AV_FILE, row.names=TRUE)
write.csv(syn, SYNERGY_FILE, row.names=FALSE)


# == DONE


# [END]