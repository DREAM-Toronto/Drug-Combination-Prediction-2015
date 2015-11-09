# XvalidatePredictions.R
#
# Purpose: From input Training Data, select a subset for
#              Training and prediction, run
#              similarity computations,
#              predict synergy scores, and evaluate
#              success.
#
# Version: 0.2.1           
#
# Date:    Nov 8 2015
# Author:  Boris and DREAM team UofT
#          
# V 0.2.1  Add analysis and plots for multiple xVal runs.
#          Handle edge cases (no data, combinations not
#              unique ... )
#          Add progress messages.
# V 0.2    Maintenance and refactoring.
#          Loop over Xval runs.
#          Calculate random correlations
# V 0.1    First code
# ==========================================================

setwd(DREAMDIR)
source("DREAMutilities.R")

# == CONSTANTS =============================================
#
METHOD <- "QuickPredict"

RUN_PATH <- "../validate/" # Directory for validation run data
RUN_PREFIX <- "08.33pct" # Filename prefix

FRAC_HOLD <- 0.33  # Fraction of holdout data

N_RUNS <- 40       # Number of cross-validation runs

EXCLUDE_POOR_QA <- TRUE  # if TRUE, ignore all experiments
                         # with QA != 1
                         
VERBOSE <- TRUE   # Print progress

# Probably no need to change below here ====================                         
# Master file that contains monotherapy data
MASTER_DATA_FILE <- "../Challenge Data/Drug Synergy Data/ch1_train_combination_and_monoTherapy.csv"


# == PACKAGES ==============================================


# == FUNCTIONS =============================================
#
xValColours <- function(x,
                        nBreaks = 10,
                        method = "highlightMaximum") {
	# create various types of coloring 
	if (method == "highlightMaximum") {
		# highlight large values with red,
		fCol <- colorRampPalette(c("#333333", "#4d4d4d",
		                           "#666666", "#808080",
		                           "#999999", "#b3b3b3",
		                           "#cccccc", "#fc9973",
		                           "#ff0d35")) 
	}
    colVec <- fCol(nBreaks)
    x <- colVec[cut(x, breaks = nBreaks, labels=FALSE)]
    return(x)	
}

# == MAIN ==================================================
#

# master file from which to select
master <- read.csv(MASTER_DATA_FILE, stringsAsFactors=FALSE)

if (EXCLUDE_POOR_QA) {
	master <- master[master[,"QA"] == 1, ]
}



# == LOOP OVER CROSS-VALIDATION RUNS =======

for (iRun in 1:N_RUNS) {
	
	# == MAKE RUN-SPECIFIC FILENAMES =======
	RUN_ID <- sprintf("%s.%02d", RUN_PREFIX, iRun)
	TRAINING_SET_FILE <- sprintf("%s%s_TrainingSet.csv", RUN_PATH, RUN_ID)
	TEST_SET_FILE <- sprintf("%s%s_HoldoutSet.csv", RUN_PATH, RUN_ID)
	PREDICTED_FILE <- sprintf("%s%s_prediction.csv", RUN_PATH, RUN_ID)
	PRIORITY_FILE <- sprintf("%s%s_combination_priority.csv", RUN_PATH, RUN_ID)

	# == MAKE TRAINING AND HOLDOUT SETS =======
	nHold <- round(nrow(master) * FRAC_HOLD) # number of holdouts
	iHold <- sample(1:nrow(master), nHold)   # random row index of holdouts
	
	trainSet <- master[-(iHold), ]
	holdSet  <- master[iHold, ]
	holdSet[, "SYNERGY_SCORE"] <- NA  # remove synergy values from
	                                  # holdout data
	
	write.csv(trainSet, TRAINING_SET_FILE, row.names=FALSE)
	write.csv(holdSet, TEST_SET_FILE, row.names=FALSE)
	
	# == RUN PREDICTION =======
	
	# predictSynergy.R is written to take input training and test
	# filenames from commandline. If it it source()'ed instead
	# redefine the commandArgs() function as below. 
	
	commandArgs <- function(trailingOnly) {
		return(c(TRAINING_SET_FILE,
		         TEST_SET_FILE,
		         PREDICTED_FILE,
		         PRIORITY_FILE,
		         METHOD))
	}
	if (VERBOSE) {print(paste("xVal run: ", iRun, " / ", N_RUNS, "    ", Sys.time()))}
	source("predictSynergy.R")
	
	# predictSynergy.R writes its result into 
	# PREDICTED_FILE and PRIORITY_FILE
	
	# == READ AND COMBINE RESULTS =======
	pred <- read.csv(PREDICTED_FILE, stringsAsFactors=FALSE)
	conf <- read.csv(PRIORITY_FILE, stringsAsFactors=FALSE)
    # head(pred)
    # head(conf)
	
	results <- data.frame(RUN = rep(iRun, nrow(pred)),
	                      CELLS  = pred[ , "CELL_LINE"],
	                      ID  = pred[ , "COMBINATION_ID"],
	                      TRUE_SYN  = NA,
	                      PRED_SYN  = pred[ , "PREDICTION"],
	                      SYN_CONF  = rep(0, nrow(pred)),
	                      stringsAsFactors = FALSE)
    for (i in 1:nrow(results)) {
    	# add confidence values to results
    	results[i, "SYN_CONF"] <- conf[ (conf[,"COMBINATION_ID"] == results[i, "ID"]) ,
    	                                "CONFIDENCE" ]
    	# add true synergy score to results
    	# use the mean, since we actually have cases where
    	# there are more than one experiments for the same
    	# combination ID and cell line.
    	results[i, "TRUE_SYN"] <- master[master[ , "COMBINATION_ID"] == pred[i, "COMBINATION_ID"] &
    	                                 master[ , "CELL_LINE"]      == pred[i, "CELL_LINE"],
    	                                     "SYNERGY_SCORE"] %>% mean
    }
    # head(results)
    
    xCor <- cor(results[,"TRUE_SYN"],
                results[,"PRED_SYN"],
                use="complete.obs")
    if (iRun == 1) {
    	allResults <- results
    	allCorrelations <- xCor
    } else {
    	allResults <- rbind(allResults, results)
     	allCorrelations <- c(allCorrelations, xCor)
   }
}   # end for (iRun in 1:N_RUNS)

if (VERBOSE) {print(paste("xVal completed: ", Sys.time()))}


# == ANALYSE RUNS =======

# == Scatterplot true vs. pred

true <- allResults[,"TRUE_SYN"]
pred <- allResults[,"PRED_SYN"]
plot(true, pred,
     cex = 0.2,
     pch = 16,
     col = xValColours(allResults[,"SYN_CONF"]),
     xlim = c(min(c(true, pred)), max(c(true, pred))),
     ylim = c(min(c(true, pred)), max(c(true, pred))),
     xlab = "True synergy scores",
     ylab = "Predicted synergy scores",
     main = sprintf("%s crossvalidation with %d%% holdouts, %d runs",
                    METHOD,
                    FRAC_HOLD * 100,
                    N_RUNS),
     cex.main = 0.8
    )
legend("topleft",
       cex = 0.6,
       legend = paste(seq(0.9, 0.0, by = -0.1), "-", seq(1.0, 0.1, by = -0.1)),
       fill = xValColours(seq(0.95, 0.05, by = -0.1)),
       bty = "n",
       title = "Priority")
    
abline(h = 0, lwd = 0.5, col = "#CCCCCC")
abline(v = 0, lwd = 0.5, col = "#CCCCCC")

# add regression line
linMod <- lm(true ~ pred)
abline(linMod, col="#BB0000")
summary(linMod)

# distribution of observed correlations
allCorrelations
mean(allCorrelations)
sd(allCorrelations)

# calculate distribution of random correlations
nRandRuns <- 1000  # per xVal run
corRand <- numeric(nRandRuns * N_RUNS)
iC <- 1
for (i in 1:N_RUNS) {
	# subset of predictions for this run
	tp <- allResults[allResults[, "RUN"] == i, c("TRUE_SYN", "PRED_SYN")]
	for (j in 1:nRandRuns) {
		# sample from subset
	    corRand[iC] <- cor(tp[ , "TRUE_SYN"],
	                      sample(tp[ , "PRED_SYN"], nrow(tp)))
	    iC <- iC + 1
	}
}


# Plot histogram of random correlations,
# superimpose histogram observed prediction correlations
colRand <- "#F5F5FF"
colObs <- "#00DDAA33"
hist(corRand, 
     breaks=20, 
     xlim = c(min(min(corRand), min(allCorrelations)) * 1.2,
              max(max(corRand), max(allCorrelations)) * 1.2),
     freq = FALSE,
     col = colRand,
     main = "Random and observed correlations",
     cex.main = 0.8,
     xlab = "R"
    )

if (length(allCorrelations) == 1) {
	abline(v = allCorrelations[rep(1,4)], col = colObs)
} else {
	par(new = TRUE)
	hVals <- hist(allCorrelations, plot = FALSE)
	hist(allCorrelations, 
	     xlim = c(min(min(corRand), min(allCorrelations)) * 1.2,
	              max(max(corRand), max(allCorrelations)) * 1.2),
	     ylim = c(0, 5 * max(hVals$density)),
	     breaks = 8,
	     freq = FALSE,
	     col= colObs,
	     border= "#666677",
	     axes = FALSE,
	     main = "",
	     xlab = "",
	     ylab = ""
	    )
}

legend("topright",
       cex = 0.7,
       legend = c(sprintf("%d,%03d random",
                          floor(length(corRand)/1000),
                          length(corRand) %% 1000),
                  sprintf("40 observed", N_RUNS)),
       fill = c(colRand, colObs),
       bty = "n",
       title = "")



# [END]
