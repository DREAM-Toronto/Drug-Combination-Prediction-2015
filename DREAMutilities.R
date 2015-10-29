# DREAMutilities.R
#
# Purpose: Utility functions for working with raw drug data
# Version: 0.3
# Date:    Oct 29 2015
# Author:  Boris and DREAM team UofT
#
# ToDo:    Make surface plot log/log concentration
#
# V 0.3    Add NLS fit for monotherapy and plot
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



# ==== nlsDRC ==============================================
# non-linear least-squares fit of a Dose Response
# Curve.

nlsDRC <- function(drc) {
    Dose   = drc$conc
    Effect = drc$dat
    nlsHill <- nls(Effect ~ fHill(Dose, IC50, H, Einf ),
                start = c(IC50 = max(Dose)/2,
                          H    = 5,
                          Einf = min(Effect)),
                lower = c(IC50 = Dose[2], H = 0,   Einf = 0),
                upper = c(IC50 = Dose[6], H = 10,  Einf = 100),
                algorithm = "port")
    return(nlsHill)
}

# ==== Hill coefficient function ===========================
# Standard 4 parameter dose-effect function with E0
# set to 100

fHill <- function(dose, IC50, H, Einf) {
	E0 <- 100
	E <- E0 + ((Einf - E0) / (1 + ((IC50 / dose) ^ H)))
    return(E)
} 


# ==== Plot DRC fit =======================================
# Calculate and plot an NLS fit for a Dose Response Curve

plotDRCfit <- function(drc) {
    fit <- nlsDRC(drc)

    scale <- drc$conc %>% reRangeConc %>% log

    plot(scale, drc$dat,
         ylim = c(0, 100),
         xlab = sprintf("log([%s])", drc$A),
         ylab = "% survivors",
         main = sprintf("Monotherapy: %s in %s", drc$A, drc$C))

    x <- seq(min(scale), max(scale), by = 0.01)
    abline(h=100-((100-coef(fit)["Einf"])/2), col="#DDDDFF")
    abline(v=log(coef(fit)["IC50"]), col="#DDDDFF")
    points(x, fHill(exp(x),
                    IC50 = coef(fit)["IC50"],
                    H =    coef(fit)["H"],
                    Einf = coef(fit)["Einf"]), 
           col="#AA0000",
           type="l")
           

    return(fit)
} 


# ==== plotDRS =============================================
# plots a DRS list in a perspective plot
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
if (FALSE) {
# skip all of this code when source'ing	
	
	
	
# Define data
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


# make prediction
predDRS <- makeAdditive(obsDRS)
synS <- calcDifference(predDRS, obsDRS)

# Compare observation and prediction
opar <- par()
par(mfrow = c(2,2)) 
plotDRS(obsDRS, col="lightblue")
plotDRS(predDRS, col="seagreen")
plotDRS(synS, col="firebrick",
        zlim=c(3/4 * min(synS$drDat), 4/3 * max(synS$drDat)))
par <- opar


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
vals <- c(scores[,"sCalc"], scores[,"sObs"])
lim <- c(min(vals), max(vals))
plot(scores, cex=0.8,
     xlim = lim, ylim = lim)
     
cor(scores[,"sCalc"], scores[,"sObs"])
reg <- lm(scores[,"sObs"] ~ scores[,"sCalc"])
abline(reg, col="firebrick")
abline(h=0, col="#DDDDFF")
abline(v=0, col="#DDDDDDFF")


} #end if (FALSE) ...
# END
