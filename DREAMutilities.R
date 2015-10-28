# DREAMutilities.R
#
# Purpose: Utility functions for working with raw drug data
# Version: 0.2.1
# Date:    Oct 22 2015
# Author:  Boris and DREAM team UofT
#
# ToDo:    Make surface plot log/log concentration
#
# V 0.2.1  Plot DRS on logarithmic scale
# V 0.2    Updated paths etc. and comitted to Repository
#          Created list object for DRS data
#          
# V 0.1    First code
# ==========================================================

# Define this variable as a path to your repository in your
# ~/.Rprofile
setwd(DREAMDIR)


# == CONSTANTS =============================================
#

# If you put your challenge data folder on the folder
# that also contains the local copy of the Repository,
# this is the right path to use. Recommended.
TRAINING <- "../Challenge Data/Drug Synergy Data/Raw Data/Raw_Data_csv/ch1_training_combinations/"


# == PACKAGES ==============================================
#
if (! require(plot3D, quietly=TRUE)) {
	# for plotting dose-response surfaces
	install.packages("plot3D",
	                 repos="http://cran.us.r-project.org")
	library(plot3D)
}

# == FUNCTIONS =============================================
#

# ==== readDRS =============================================
readDRS <- function(fn, path=TRAINING) {
	# Read raw drug interaction data from csv file containing
	# a dose-response surface (DRS).
	# Arguments
	#   fn: filename
	#   path: default TRAINING.
	# We assume row 1, column 1 are the concentrations
	# and the data is in a 6x6, fully filled matrix
	# of percent-surviving cells.
	# Return a list that contains:
	#    $A:     name of compound A
	#    $B:     name of compound B
	#    $C:     name of cell line
	#    $concA: concentrations for compound A
	#    $concB: concentrations for compound B
	#    $drDat: dose response data (%survivors)
	drs <- list()
	raw <- read.csv(paste(path, fn, sep=""),
                    head=FALSE,
                    stringsAsFactors=FALSE)
    drs$A <- raw[ 9, 2]
    drs$B <- raw[10, 2]
    drs$C <- raw[13, 2]
    drs$type <- "% survivors"
    drs$concA <- as.numeric(raw[2:7, 1])
    drs$concB <- as.numeric(raw[1, 2:7])
    drs$drDat <- raw[2:7, 2:7] %>%
           unlist %>%
           as.numeric %>%
           as.matrix
    dim(drs$drDat) <- c(6, 6)
    rownames(drs$drDat) <- drs$concA
    colnames(drs$drDat) <- drs$concB
    return(drs)
}


# ==== makeFileName ========================================
# Make a file-name string from compound A, compound B
# and cell line name.

makeFileName <- function(a, b, c) {
    return(paste(a, b, c, "Rep1.csv", sep="." ))
}


# ==== getMono =============================================
# Get a monotherapy DRC list from a DRS.

getMono <- function(DRS, compound) {
	DRC <- list()
	DRC$A <- compound
	DRC$C <- DRS$C
	if (DRC$A == DRS$A) {
		DRC$conc <- DRS$concA
		DRC$dat  <- DRS$drDat[, "0"]
	} else if (DRC$A == DRS$B) {
		DRC$conc <- DRS$concB
		DRC$dat  <- DRS$drDat["0", ]
	} else {
		stop("Requested compound not in DRS list.")
	}
    return(DRC)
}


# ==== makeAdditive ========================================
makeAdditive <- function(drs) {
# Makes a combination response matrix with a simple
# additive model.
# Returns a DRS list
	add <- drs
	add$type <- "% survivors (additive prediction)"
	for (i in 1:6) {
		for (j in 1:6) {
			addVal <- 100 - ((100-drs$drDat[i,1]) + (100-drs$drDat[1,j]))
			add$drDat[i, j] <- max(addVal, 0)  # set negative values to 0
		}
	}
	return(add)
}

# ==== calcDifference ======================================
calcDifference <- function(a, b) {
# Calculate the difference between data matrices in two
# DRS lists. This is probably only reasonable for matrices
# that represent the same compounds and cells, such as
# comparing the predicted additive vs. the observed effect
# on survival, to establish synergy; however this function
# does not enforce this, nor does it enforce that the 
# concentration ranges are identical.
# Return a DRS list.
    diff <- a
	diff$type <- "% DRSa - DRSb"
	diff$drDat <- a$drDat - b$drDat
	return(diff)
}


# ==== integrateDRS ========================================
# Integrates a DRS matrix
# Right now this is simply the sum, but we might do
# something more sophisticated in the future ...

integrateDRS <- function(drs) {
    return(sum(drs$drDat))
}


# ==== plotDRS =============================================
# plots a DRS list
plotDRS <- function(drs, col="white", zlim= c(0,110)) {
  scaleA <- drs$concA %>% reRangeConc %>% log
  scaleB <- drs$concB %>% reRangeConc %>% log
  persp3D(x = scaleA,
        y = scaleB,
        z = drs$drDat,
        theta = 135, phi = 30,
        xlim=c(0.9*scaleA[1], 1.1*scaleA[6]),
        ylim=c(0.9*scaleB[1], 1.1*scaleB[6]),
        zlim=zlim,
        col = col,
        border = "black",
        ltheta = 120, shade = 0.75,
        ticktype = "detailed",
        cex.axis = 0.67,
        cex.lab = 0.67,
        cex.main = 0.80,
        xlab = paste("log [", drs$A, "]", sep=""),
        ylab = paste("log [", drs$B, "]", sep=""),
        zlab = drs$type,
        main = paste(drs$A, "+", drs$B, "in", drs$C))
}

reRangeConc <- function(v) {
	# re-range a 6-element vector of compound
	# concentrations to set the first element
	# (which is 0) to one log step smaller than
	# the second element, so the data can be 
	# displayed on a log scale.
	fac <- v[3] / v[2]
	v[1] <- v[2] / fac
	return(v)
}


# ==== explore ========================================

compoundA <- "ADAM17"
compoundB <- "AKT"
cells <- "HCC1806"
fn <- makeFileName(compoundA, compoundB, cells)

obsDRS <- readDRS(fn)
predDRS <- makeAdditive(obsDRS)
synS <- calcDifference(predDRS, obsDRS)

plotDRS(obsDRS, col="lightblue")
plotDRS(predDRS, col="seagreen")
plotDRS(synS, col="firebrick",
        zlim=c(3/4 * min(synS$drDat), 4/3 * max(synS$drDat)))


iObs <- integrateDRS(obsDRS)
iPred <- integrateDRS(predDRS)
iSyn <- integrateDRS(synS)

synScore <- 100 * (iPred - iObs) / iPred
synScore

# Note that this value does not implement the _actual_
# Lowe model, and integration is not done in log-dose
# space. This needs to be refined.

# Test our simple synScores against the values found in
# ch1_train_combination_and_monoTherapy.csv

fname <- "../Challenge Data/Drug Synergy Data/ch1_train_combination_and_monoTherapy.csv"

train <- read.csv(fname, stringsAsFactors=FALSE)
head(train)
nrow(train)

# compute pred. and obs for n random rows
n <- 300
set.seed(112358)
randRows <- sample(nrow(train), n)
scores <- matrix(numeric(2*n), nrow=n, ncol=2)
colnames(scores) <- c("sCalc", "sObs")
ind <- 0
for (i in 1:n) {
	row <- randRows[i]
	if (train[row,"QA"] == 1) {
		A <- train[row,"COMPOUND_A"]
		B <- train[row,"COMPOUND_B"]
		C <- train[row,"CELL_LINE"]
	    fn <- makeFileName(A, B, C)
	    obsDRS <- readDRS(fn)
		predDRS <- makeAdditive(obsDRS)
		synS <- calcDifference(predDRS, obsDRS)
        iObs <- integrateDRS(obsDRS)
        iPred <- integrateDRS(predDRS)
        ind <- ind + 1
        scores[ind, "sCalc"] <- 100 * (iPred - iObs) / iPred
        scores[ind, "sObs"] <- train[row,"SYNERGY_SCORE"]
	}
}

scores <- scores[1:ind, ]
plot(scores, cex=0.8)

corr <- lm(scores[,"sCalc"] ~ scores[,"sObs"])
abline(corr, col="firebrick")

# END
