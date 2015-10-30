# DREAMexploreData.R
#
# Purpose: Sample explorations of raw drug data
# Version: 0.1
# Date:    Oct 30 2015
# Author:  Boris and DREAM team UofT
#
# V 0.1    Code separated out from DREAMutilities.R V0.4.1
#
# ==========================================================

# Variables and functions used here require that
# DREAMutilities.R should have been loaded. If this is setup
# as recommended, the file has been sourced()'d by your
# ~/.Rprofile stratup script.

setwd(DREAMDIR)

# Look at a single DRS (Dose Response Surface)	
# Define compounds and cells
compoundA <- "ADAM17"
compoundB <- "MTOR_1"
cells <- "DU-4475"
fn <- makeFileName(compoundA, compoundB, cells)

# read data
obsDRS <- readDRS(fn)

# plot monotherapy
drcA<- getMono(obsDRS, "ADAM17")
plotDRCfit(drcA)

drcB<- getMono(obsDRS, "MTOR_1")
plotDRCfit(drcB)


# plot Dose Response Surface
plotDRS(obsDRS)


# make sequential additive prediction
predDRS <- makeAdditive(obsDRS)
synS <- calcDifference(predDRS, obsDRS)

# Compare observation, additive effect, 
# and "synergy" in the same window
opar <- par(mfrow = c(2,2), omi = c(0, 0, 0, 0), mar = c(1, 0.5, 1, 0.5)) 
plotDRS(obsDRS, col="lightblue")
plotDRS(predDRS, col="seagreen")
plotDRS(synS, col="firebrick",
        zlim=c(3/4 * min(synS$drDat), 4/3 * max(synS$drDat)))
par <- opar


# Calculate a synergy score by integrating in log-space
iObs <- integrateDRS(obsDRS)
iObs
iPred <- integrateDRS(predDRS)
iPred
iSyn <- integrateDRS(synS)
iSyn

synScore <- 100 * (iPred - iObs) / iPred
synScore


# Is our calculated synergy score the same that the 
# Challenge team reports? Compare our synScores against
# the values reported in ch1_train_combination_and_monoTherapy.csv

fname <- "../Challenge Data/Drug Synergy Data/ch1_train_combination_and_monoTherapy.csv"

# read the training data
train <- read.csv(fname, stringsAsFactors=FALSE)
head(train)
nrow(train)


# compute pred. and reported for n random rows
n <- 200
set.seed(123456)

randRows <- sample(nrow(train), n) # pick n rows at random

scores <- matrix(numeric(2*n), nrow=n, ncol=2)
colnames(scores) <- c("synCalc", "synRep")

ind <- 0  # index of rows for which calculations were done 
for (i in 1:n) {
	row <- randRows[i]
	if (train[row,"QA"] == 1) { # make a prediction
		A <- train[row,"COMPOUND_A"]
		B <- train[row,"COMPOUND_B"]
		C <- train[row,"CELL_LINE"]
	    fn <- makeFileName(A, B, C)
	    obsDRS <- readDRS(fn)  # actual observed DRS
		predDRS <- makeAdditive(obsDRS) # additive prediction
        iObs <- integrateDRS(obsDRS)
        iPred <- integrateDRS(predDRS)
        ind <- ind + 1
        scores[ind, "synCalc"] <- 100 * (iPred - iObs) / iPred
        scores[ind, "synRep"] <- train[row,"SYNERGY_SCORE"]
	}
}

# drop empty lines (QA was not 1)
scores <- scores[1:ind, ]

# plot
plot(scores, cex=0.8, main = "Reported vs. calculated synergy")

# calculate correlation and plot regression line
cor(scores[,"synCalc"], scores[,"synRep"])
reg <- lm(scores[,"synRep"] ~ scores[,"synCalc"])
abline(reg, col="#CC0000")
abline(h=0, col="#DDDDFF")
abline(v=0, col="#DDDDFF")

# Bottom line: our approximation and the values calculated
# by Combenefit are similar (up to scale), the differences
# don't seem to be systematic but due to noise in the
# calculation.

# Let's do the same calculation for experiments for which
# both compounds had an Einf < 70 and an IC50 within the
# concentration range...

n <- 500
set.seed(123456)

randRows <- sample(nrow(train), n) # pick n rows at random

scores <- matrix(numeric(2*n), nrow=n, ncol=2)
colnames(scores) <- c("synCalc", "synRep")

ind <- 0  # index of rows for which calculations were done 
for (i in 1:n) {
	row <- randRows[i]
	    
	if (train[row,"QA"] == 1 &&
	    train[row,"IC50_A"] < train[row,"MAX_CONC_A"] &&
	    train[row,"IC50_A"] > train[row,"MAX_CONC_A"]/100 &&
	    train[row,"Einf_A"] < 70 &&
	    train[row,"IC50_B"] < train[row,"MAX_CONC_B"] &&
	    train[row,"IC50_B"] > train[row,"MAX_CONC_B"]/100 &&
	    train[row,"Einf_B"] < 70
	    ) {
		A <- train[row,"COMPOUND_A"]
		B <- train[row,"COMPOUND_B"]
		C <- train[row,"CELL_LINE"]
	    fn <- makeFileName(A, B, C)
	    obsDRS <- readDRS(fn)  # actual observed DRS
		predDRS <- makeAdditive(obsDRS) # additive prediction
        iObs <- integrateDRS(obsDRS)
        iPred <- integrateDRS(predDRS)
        ind <- ind + 1
        scores[ind, "synCalc"] <- 100 * (iPred - iObs) / iPred
        scores[ind, "synRep"] <- train[row,"SYNERGY_SCORE"]   	
	}
}

ind   # How many did we find?


# drop empty lines (QA was not 1)
scores <- scores[1:ind, ]

# plot
plot(scores, cex=0.8, main = "Reported vs. calculated synergy")

# calculate correlation and plot regression line
cor(scores[,"synCalc"], scores[,"synRep"])
reg <- lm(scores[,"synRep"] ~ scores[,"synCalc"])
abline(reg, col="#CC0000")
abline(h=0, col="#DDDDFF")
abline(v=0, col="#DDDDFF")

# Result: the noise in the comparison is not due to the
# experiments being "better" or "worse".




# END
