# makeFeatureVectors.R
# Produces feature vector and synergy data files from
# training data.
#
# Version: 0.4           
#
# Date:    Nov 13 2015
# Author:  Boris and DREAM Toronto
#
# V 0.4    Fixes incorrect handling of normalization by
#              normalizing drug IC50 to max(MAX_CONC)
#          Add PCA reduction of feature numbers
#          Add a bit of analysis, plotting
# V 0.3    Add output files for full-length combination
#              feature vectors.
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
OUT_COMBI_FILE <- "combiFeaturesAv.csv" # concatenated features for combinations
OUT_COMBI_NO_H_FILE <- "combiFeaturesAvNoH.csv" # as above but no H features
OUT_COMBI_PC_FILE <- "combiFeaturesPC.csv"  # most predictive Principal Components only
OUT_SYN_FILE <- "drugSynergies.csv"

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


# == MAX_CONC VALUES TO NORMALIZE DRC PARAMETERS

# matrix of concentration data:
# nExp <- 300 is enough for the ch1 training set - if more
# experiments for one drug are ever encountered, this part
# will fail with a subscript out of bounds error. Then
# increase the parameter.
nExp <- 300
conc <- matrix(numeric(length(drugs) * nExp), ncol=nExp)
conc[, 1] <- 2 # put an index for next empty slot of each row
               # into column 1. We'll update the index as we
               # find the various experiments' MAX_CONC
               # values in "trainingData".
rownames(conc) <- drugs

for (i in 1:nrow(trainingData)) {
    A <- trainingData$COMPOUND_A[i]
    B <- trainingData$COMPOUND_B[i]
	conc[A, conc[A, 1]] <- trainingData$MAX_CONC_A[i]
	conc[A, 1] <- conc[A, 1] + 1  # increment slot index
	conc[B, conc[B, 1]] <- trainingData$MAX_CONC_B[i]
	conc[B, 1] <- conc[B, 1] + 1
}
# summary(conc[ , 1]) # numbers of experiments per drug
# head(conc)

# tabulate results
# varConc <- 0
# for (i in 1:nrow(conc)) {
	# if (length(unique(conc[i, 2:(conc[i, 1] - 1)])) > 1) {
		# print(paste(rownames(conc)[i], paste(conc[i, 2:(conc[i, 1] - 1)], collapse=" ")))
		# varConc <- varConc + 1
	# }
# }
# print(sprintf("Differing concentrations in %3.1f%% of drugs", 100 * varConc / nrow(conc)))

# get the row maxima from the matrix
maxConc <- apply(conc[ , 2:nExp], MARGIN = 1, FUN = max)
# head(maxConc)

# normalize concentrations and H
for (i in 1:nrow(trainingData)) {
    A <- trainingData $COMPOUND_A[i]
    B <- trainingData $COMPOUND_B[i]

    trainingData$IC50_A[i] <- trainingData$IC50_A[i] / maxConc[A]
    trainingData$IC50_B[i] <- trainingData$IC50_B[i] / maxConc[B]

    trainingData$H_A[i] <- trainingData$H_A[i] * maxConc[A]
    trainingData$H_B[i] <- trainingData$H_B[i] * maxConc[B]

}
# head(trainingData)


# == PREPARE DRUG AND CELL LINE FEATURE MATRIX

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


# Make a copy of data in which we replace NA values 
# with row averages. This is probably not the best way
# to impute these missing values. ToDo: Look into the
# imputation literature...

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

# make the combi-features datasets
combiFeatures <- matrix(numeric(nrow(syn) * ncol(data) * 2), 
                        nrow = nrow(syn))
rownames(combiFeatures) <- syn[ , "COMBINATION_ID"]
colnames(combiFeatures) <- rep(colnames(data), 2)
                        

iTriplet <- seq(1, ncol(combiFeatures), by = 3)

for (i in 1:nrow(syn)) {
	# concatenate compound A and B features
	combiFeatures[i, ] <- c(data2[syn[i, "COMPOUND_A"], ],
	                        data2[syn[i, "COMPOUND_B"], ])
	# scale features to zero mean and unit variance
	combiFeatures[i, iTriplet]   <- scale(combiFeatures[i, iTriplet])
	combiFeatures[i, iTriplet+1] <- scale(combiFeatures[i, iTriplet+1])
	combiFeatures[i, iTriplet+2] <- scale(combiFeatures[i, iTriplet+2])
}
# head(combiFeatures)

# make version without H values
iH <- seq(2, ncol(combiFeatures), by = 3)
combiFeaturesNoH <- combiFeatures[ , -iH]
# head(combiFeaturesNoH)


# == DIMENSION REDUCTION

# Part 1: which of the triplet features seem important
#         for the correlation of synergy scores with features?
#         We correlate the correlation of drug-pair features 
#         with the drug-pair synergy scores. This is a 
#         correlation of correlations.

datNew <- data2

# make vectors of column indices for the DRC parameters
iIC <- seq(1, ncol(datNew), by = 3)
iH  <- seq(2, ncol(datNew), by = 3)
iEi <- seq(3, ncol(datNew), by = 3)

# scale columns to 0 mean and unit variance
for (i in 1:ncol(datNew)) {
	datNew[ , i] <- scale(datNew[ , i])
}

corNew <- matrix(numeric(nrow(syn) * 7), ncol = 7)
# calculate drug-pair correlations for each feature
# correlation
colnames(corNew) <- c("IC", "H", "Ei", "IC.H", "IC.Ei", "H.Ei", "IC.H.Ei")

selCor <- function(A, B, columns) {
	# calculate row correlations for drugs A, B in selected columns only
	return(cor(datNew[A, columns], datNew[B, columns]))
}


for (i in 1:nrow(syn)) {
	A <- syn[i, "COMPOUND_A"]
	B <- syn[i, "COMPOUND_B"]
	corNew[i, "IC"     ] <-  selCor(A, B, iIC)
	corNew[i, "H"      ] <-  selCor(A, B, iH)
	corNew[i, "Ei"     ] <-  selCor(A, B, iEi)
	corNew[i, "IC.H"   ] <-  selCor(A, B, c(iIC, iH ))
	corNew[i, "IC.Ei"  ] <-  selCor(A, B, c(iIC, iEi))
	corNew[i, "H.Ei"   ] <-  selCor(A, B, c(iH,  iEi))
	corNew[i, "IC.H.Ei"] <-  selCor(A, B, c(iIC, iH,  iEi))
}

cor(syn$SYNERGY_SCORE, corNew[ , "IC"     ])  # -0.0220
cor(syn$SYNERGY_SCORE, corNew[ , "H"      ])  #  0.0696
cor(syn$SYNERGY_SCORE, corNew[ , "Ei"     ])  # -0.0897
cor(syn$SYNERGY_SCORE, corNew[ , "IC.H"   ])  #  0.1360
cor(syn$SYNERGY_SCORE, corNew[ , "IC.Ei"  ])  #  0.0217
cor(syn$SYNERGY_SCORE, corNew[ , "H.Ei"   ])  #  0.1489  ***
cor(syn$SYNERGY_SCORE, corNew[ , "IC.H.Ei"])  #  0.1360

# This shows that all three parameters contain useful
# information ONLY in combination with each other.
# Even though "H.Ei" is a bit better than "IC.H.Ei"
# we'll continue with all three parameters and let PCA
# sort out how best to combine them.

plot(syn$SYNERGY_SCORE,
     corNew[ , "IC.H.Ei"],
     main="Mean synergy score vs. drug combination parameter correlation",
     xlab="mean Synergy Score for A.B",
     ylab="cor(A, B)",
     cex=0.8,
     cex.main=0.7,
     cex.lab=0.8,
     pch=21,
     bg="#DDEEFF")
mod <- lm(corNew[ , "IC.H.Ei"] ~ syn$SYNERGY_SCORE)
abline(mod, col="#BB0000")
text(-20,0.7, "R = 0.136", cex=0.8)



# Part 2: Calculate PCA and select PCs

pcaDrugs <- prcomp(datNew, na.rm=TRUE)
# plot(pcaDrugs)
# summary(pcaDrugs) # ~75% of variance in the first five PCs
# str(pcaDrugs)


# Which number of PCs gives the best correlation with synergy score?
x <- numeric()
for (ii in 2:ncol(pcaDrugs$x)) { # for PCs 1 to ii
    corPC <- numeric(nrow(syn))
    iCol <- 1:ii
    for (i in 1:nrow(syn)) {
        corPC[i] <- cor(pcaDrugs$x[syn[i, "COMPOUND_A"], iCol],
	                    pcaDrugs$x[syn[i, "COMPOUND_B"], iCol]) 
    }
    x[ii] <- cor(syn$SYNERGY_SCORE, corPC)
}
plot(x, ylim = c(0, 0.12))  # best choice at ii == 5

nBest <- 5

# Plot the correlation
corPC <- numeric(nrow(syn))
iCol <- 1:nBest
for (i in 1:nrow(syn)) {
    corPC[i] <- cor(pcaDrugs$x[syn[i, "COMPOUND_A"], iCol],
                    pcaDrugs$x[syn[i, "COMPOUND_B"], iCol]) 
}
cor(syn$SYNERGY_SCORE, corPC)

plot(syn$SYNERGY_SCORE,
     corPC, 
     main="Mean synergy score vs. drug combination parameter\n correlation of first five Principal Components",
     xlab="mean Synergy Score for A.B",
     ylab="cor(PC$x(A, B)[1:nBest])",
     cex=0.8,
     cex.main=0.7,
     cex.lab=0.8,
     pch=21,
     bg="#DDEEFF")
mod <- lm(corPC ~ syn$SYNERGY_SCORE)
abline(mod, col="#BB0000")
text(-18, 0.85, "R = 0.100", cex=0.8)

cRand <- numeric(10000)
for (i in 1:10000) {
	cRand[i] <- cor(syn$SYNERGY_SCORE, corPC[sample(1:length(corPC))])
}

hist(cRand, breaks=25, col="#F7FDFE", main="Mean synergy score and drug combination parameter\n correlation of first five PCs: random vs. predicted", xlab="R", cex.main=0.7)
abline(v=cor(syn$SYNERGY_SCORE, corPC), col="#BB0000", lwd=1.5)
# Poor significance. But that's what it is.



# Part 3: turn nBest PCs into a combi-features dataset
combiFeaturesPC <- matrix(numeric(nrow(syn) * nBest * 2), 
                          nrow = nrow(syn))
rownames(combiFeaturesPC) <- syn[ , "COMBINATION_ID"]
colnames(combiFeaturesPC) <- c(paste("PCA", 1:nBest, sep="."), 
                               paste("PCB", 1:nBest, sep=".") )
                        
iCol <- 1:nBest
for (i in 1:nrow(syn)) {
	# concatenate compound A and B features
	combiFeaturesPC[i, ] <- c(pcaDrugs$x[syn[i, "COMPOUND_A"], iCol],
                              pcaDrugs$x[syn[i, "COMPOUND_B"], iCol]) 
}
# head(combiFeaturesPC)


# == WRITE FILES
write.csv(data,  OUT_NA_FILE, row.names=TRUE)
write.csv(data2, OUT_AV_FILE, row.names=TRUE)
write.csv(syn, OUT_SYN_FILE, row.names=FALSE)
write.csv(combiFeatures, OUT_COMBI_FILE, row.names=TRUE)
write.csv(combiFeaturesNoH, OUT_COMBI_NO_H_FILE, row.names=TRUE)
write.csv(combiFeaturesPC, OUT_COMBI_PC_FILE, row.names=TRUE)


# == DONE


# [END]