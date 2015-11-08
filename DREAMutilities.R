# DREAMutilities.R
#
# Purpose: Utility functions for working with raw drug data
# Version: 0.4.2
# Date:    Oct 30 2015
# Author:  Boris and DREAM team UofT
#
# V 0.4.2  Changed path of training data directory to
#              create itself from DREAMDIR, defined in .Rprofile
# V 0.4.1  Improve logic of nlsDRC, make fit more robust.
#          Separate out data exploration examples into
#              DREAMexploreData.R
# V 0.4    Handle convergence failures in NLS fit.
#          Implement sequential additive model
#          to estimate additive compound effects.
# V 0.3    Add NLS fit for monotherapy and plot
# V 0.2.1  Plot DRS on logarithmic scale
# V 0.2    Updated paths etc. and comitted to Repository
#          Created list object for DRS data
#          
# V 0.1    First code
# ==========================================================

# Define DREAMDIR as a path to your project directory
# in your ~/.Rprofile. DREAMDIR should contain:
# - the folder "Drug-Combination-Prediction-2015", which is
#       the repository on github
# - the folder "Challenge Data", which contains 
#       -- "Drug Synergy Data" and 
#       -- "Sanger Molecular Data"
# - other folders as needed.

setwd(DREAMDIR)


# == CONSTANTS =============================================
#

TRAIN_DIR <- paste(DREAMDIR,
                   "Challenge Data/Drug Synergy Data/Raw Data/Raw_Data_csv/ch1_training_combinations/",
                   sep="")


# == PACKAGES ==============================================
#
if (! require(plot3D, quietly=TRUE)) {
	# for plotting dose-response surfaces
	install.packages("plot3D",
	                 repos="http://cran.us.r-project.org")
	library(plot3D)
}

if (! require(magrittr, quietly=TRUE)) {
    # for plotting dose-response surfaces
    install.packages("magrittr",
                     repos="http://cran.us.r-project.org")
    library(magrittr)
}


# == FUNCTIONS =============================================
#

# ==== readDRS =============================================
readDRS <- function(fn, path = TRAIN_DIR) {
	# Read raw drug interaction data from csv file containing
	# a dose-response surface (DRS).
	# Arguments
	#   fn: filename
	#   path: default TRAIN_DIR but can be defined differently
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
# Sequential additive model for drug combinations.
#
# The mathematics of Loewe's dose-additive model
# is not applicable to our situation of partial
# agonists, because an "equivalent dose" Beq of a
# compound B that has a lower Einf than a compound
# A is not defined for combinations for which
# the mono-effect of A is greater than Einf of B.
#
# I use a sequential additive model instead: 
# A compound B exerts it's effect in addition to A
# by reducing the survivors of A by the percentage
# expected from B.
#
# For a dosage pair dA, dB calculate the relative
# effect EdB, EdB individually. The combined effect is
# EdA * EdB/100 (with effects in % surviving cells).
#
# Returns a DRS list
	add <- drs
	add$type <- "% survivors (sequential additive model)"

    drcA <- getMono(drs, drs$A)
    drcB <- getMono(drs, drs$B)

    coefA <- nlsDRC(drcA)
    coefB <- nlsDRC(drcB)

	for (i in 1:6) {
		add$drDat[i, 1] <- fHill(dose = drcA$conc[i],
		                         IC50 = coefA["IC50"],
		                         H    = coefA["H"],
		                         Einf = coefA["Einf"])
		add$drDat[1, i] <- fHill(dose = drcB$conc[i],
		                         IC50 = coefB["IC50"],
		                         H    = coefB["H"],
		                         Einf = coefB["Einf"])
    }

	for (i in 2:6) {
		for (j in 2:6) {
			addVal <- add$drDat[i,1] * (add$drDat[1,j]/100)
			add$drDat[i, j] <- addVal
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
# Integrates a DRS matrix in log concentration
# space. We are simply multiplying the dx and
# dy values from each datapoint midway to
# its neighbor, by the value of the point.
# Thus this is a simple stepwise approximation,
# no interpolation is done.
integrateDRS <- function(drs) {

	dRow <- makeDelta(conc2log(drs$concA))
	dCol <- makeDelta(conc2log(drs$concB))

    sum <- 0
    for (i in 1:5) {
    	for (j in 1:5) {    		
    		sum <- sum + (dRow[i] * dCol[j] * drs$drDat[i+1, j+1])
    	}
    }
    return(sum)
}

# Utility function for integrate DRS
makeDelta <- function(v) {
	# calculates interval represented
	# by each concentration value in v
	dv <- numeric(5)
	dv[1] <- (v[2] - v[1]) / 2
	dv[2] <- (v[3] - v[1]) / 2
	dv[3] <- (v[4] - v[2]) / 2
	dv[4] <- (v[5] - v[3]) / 2
	dv[5] <- (v[5] - v[4]) / 2
    return(dv)
}


# Utility function for integrate DRS
conc2log <- function(conc) {
	# Rescaling concentrations.
	# This assumes that the concentration spans
	# two orders of magnitude and the non-zero
	# values are between 2 and 6 in a six-element
	# vector. Min and Max are set to 1 and 100,
	# and the intermediate concentrations are
	# rescaled into this range.

	conc <- conc[-1]  # drop 0 value
	scale <- (100 - 1) / (conc[5] - conc[1])
	conc <- ((conc - conc[1]) * scale) + 1 
	return(log(conc))
}


# ==== nlsDRC ==============================================
# non-linear least-squares fit of a Dose Response
# Curve. Return the coefficients of the fit, unless
# requesting the full fit. If a three-parameter fit
# is unsuccessful, try fitting only H and Einf,
# or only Einf. Note that effects are being constrained
# between 0 and 100.

nlsDRC <- function(drc, coefOnly=TRUE) {
    Dose   = drc$conc
    Effect = drc$dat
    
    # Constrain effects between 0 and 100
    Effect[Effect <   0] <- 0
    Effect[Effect > 100] <- 100

    # try three-parameter fit
    fit <- nlsDRC.3(Dose, Effect)
    
    if (fit$failed) {
    	# try two-parameter fit with IC50 <- Dose[6]
        fit <- nlsDRC.2(Dose, Effect, Dose[6])
    }

    if (fit$failed) {
    	# try two-parameter fit with IC50 <- Dose[1]
        fit <- nlsDRC.2(Dose, Effect, Dose[1])
    }

    if (fit$failed) {
    	# try one-parameter fit with IC50 <- Dose[6]
    	# (pseudo - linear)
        fit <- nlsDRC.1(Dose, Effect, IC50 = Dose[6], H = 20)
    }

    if (fit$failed) {
    	# try one-parameter fit with IC50 <- Dose[6]
    	# (super-effective)
        fit <- nlsDRC.1(Dose, Effect, IC50 = reRangeConc(Dose)[1], H = 0)
    }

    if (fit$failed) {
    	# Still failed? Give up.
    	browser()
        stop("NLS fit did not converge. Check the data.")
    }

    if (coefOnly) {
    	return(fit$coef)
    } else {
        return(fit$result)
    }	
}


# Utility function for nlsDRC
nlsDRC.3 <- function(Dose, Effect) {
	# three-parameter fit of Hill curve
	fit <- list()
    try(fit$result <- nls(Effect ~ fHill(Dose, IC50, H, Einf ),
                start = c(IC50 = max(Dose)/2,
                          H    = 5,
                          Einf = min(Effect)),
                lower = c(IC50 = Dose[2], H = 0,   Einf = 0),
                upper = c(IC50 = Dose[6], H = 10,  Einf = 100),
                algorithm = "port"),
         silent = TRUE)

    if (length(fit$result) == 0) {
    	fit$failed <- TRUE
    } else {
    	fit$failed <- FALSE
    	fit$coef <- c(coef(fit$result)["IC50"],
                      coef(fit$result)["H"],
                      coef(fit$result)["Einf"])   	
    }
    return(fit)
}


# Utility function for nlsDRC
nlsDRC.2 <- function(Dose, Effect, IC50) {
	# two-parameter fit of Hill curve: IC50 fixed; fit H, Einf
	fit <- list()
    try(fit$result <- nls(Effect ~ fHill(Dose, IC50, H, Einf ),
                start = c(H    = 5,
                          Einf = min(Effect)),
                lower = c(H = 0,   Einf = 0),
                upper = c(H = 10,  Einf = 100),
                algorithm = "port"),
         silent = TRUE)

    if (length(fit$result) == 0) {
    	fit$failed <- TRUE
    } else {
    	fit$failed <- FALSE    	
    	fit$coef <- c(IC50 = IC50,
                      coef(fit$result)["H"],
                      coef(fit$result)["Einf"])   	
    }
    return(fit)
}


# Utility function for nlsDRC
nlsDRC.1 <- function(Dose, Effect, IC50, H) {
	# one parameter fit of Hill curve: fit only Einf
	fit <- list()
    try(fit$result <- nls(Effect ~ fHill(Dose, IC50, H, Einf ),
                start = c(Einf = min(Effect)),
                lower = c(Einf = 0),
                upper = c(Einf = 100),
                algorithm = "port"),
         silent = TRUE)

    if (length(fit$result) == 0) {
    	fit$failed <- TRUE
    } else {
    	fit$failed <- FALSE    	
    	fit$coef <- c(IC50 = IC50,
                      H =    H,
                      coef(fit$result)["Einf"])   	
    }
    return(fit)
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

    coefFit <- nlsDRC(drc)

    scale <- drc$conc %>% reRangeConc %>% log

    plot(scale, drc$dat,
         ylim = c(-10, 110),
         xlab = sprintf("log([%s])", drc$A),
         ylab = "% survivors",
         main = sprintf("Monotherapy: %s in %s", drc$A, drc$C))

    x <- seq(min(scale), max(scale), by = 0.01)
    abline(h=100-((100-coefFit["Einf"])/2), col="#DDDDFF")
    abline(v=log(coefFit["IC50"]), col="#DDDDFF")
    points(x, fHill(exp(x),
                    IC50 = coefFit["IC50"],
                    H =    coefFit["H"],
                    Einf = coefFit["Einf"]), 
           col="#AA0000",
           type="l")
                 
    return(coefFit)
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

# END
