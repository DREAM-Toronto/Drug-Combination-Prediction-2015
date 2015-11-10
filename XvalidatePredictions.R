# XvalidatePredictions.R
#
# Purpose: From input Training Data, select a subset for
#              Training and prediction, run
#              similarity computations,
#              predict synergy scores, and evaluate
#              success.
#
# Version: 0.3           
#
# Date:    Nov 9 2015
# Author:  Boris and DREAM team UofT
#          
# V 0.3    Add scoring according to AZ Challenge organizers -
#              with partial correlations to remove drug
#              combination and cell-line medians.
#          Reorganize results
# V 0.2.2  Maintenance
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
source("ch1scoring_functions.R")

# == CONSTANTS =============================================
#
METHOD <- "QuickPredict"

RUN_PATH <- "../validate/" # Directory for validation run data
RUN_PREFIX <- "09.test" # Filename prefix

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

# == Datastructures to store results
combVec = c(meanR   = numeric(),
            steR    = numeric(),
            meanMSE = numeric(),
            steMSE  = numeric(),
            N       = numeric())
scores <- list(RUN = numeric(),
               COR_PLAIN  = numeric(),
               GLOBAL_SCORE  = numeric(),
               COMB_SCORE_ALL = combVec,
               COMB_SCORE_30  = combVec,
               COMB_SCORE_20  = combVec,
               COMB_SCORE_10  = combVec)

allResults <- data.frame(RUN = numeric(),
                         CELLS  = character(),
                         ID  = character(),
                         TRUE_SYN  = numeric(),
                         PRED_SYN  = numeric(),
                         SYN_CONF  = numeric(),
                         stringsAsFactors = FALSE)


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
	if (VERBOSE) {cat(paste("xVal run: ", iRun, "/", N_RUNS, "    ", Sys.time(), "\n"))}


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
    	allResults <- rbind(allResults, results)

    # head(results)
    
    # == SCORE THE RESULTS  =======
	scores$RUN <- c(scores$RUN, iRun)
    scores$COR_PLAIN  <- c(scores$COR_PLAIN,
                           cor(results[,"TRUE_SYN"],
                               results[,"PRED_SYN"],
                               use="complete.obs"))
    scores$GLOBAL_SCORE  <- c(scores$GLOBAL_SCORE,
                              getGlobalScore_ch1(results))
    scores$COMB_SCORE_ALL = rbind(scores$COMB_SCORE_ALL,
                                  getDrugCombiScore_ch1(results, conf, topX=100))
    scores$COMB_SCORE_30  = rbind(scores$COMB_SCORE_30,
                                  getDrugCombiScore_ch1(results, conf, topX=30))
    scores$COMB_SCORE_20  = rbind(scores$COMB_SCORE_20,
                                  getDrugCombiScore_ch1(results, conf, topX=20))
    scores$COMB_SCORE_10  = rbind(scores$COMB_SCORE_10,
                                  getDrugCombiScore_ch1(results, conf, topX=10))
   
	                      
   
}   # end for (iRun in 1:N_RUNS)

if (VERBOSE) {cat(paste("xVal completed: ", Sys.time(), "\n\n"))}


# == ANALYSE RUNS =======

# == HISTOGRAM RANDOM VS. PREDICTED
# create distribution of random scores
nRandRuns <- 1000  # per xVal run
scoreRand <- numeric(nRandRuns * N_RUNS)
iC <- 1
for (i in 1:N_RUNS) {
	# subset of predictions for this run
	tp <- allResults[allResults[, "RUN"] == i, c("CELLS", "ID", "TRUE_SYN", "PRED_SYN")]
	for (j in 1:nRandRuns) {
		# sample from subset
		tp[ , "PRED_SYN"] <- sample(tp[ , "PRED_SYN"], nrow(tp))
	    scoreRand[iC] <- getGlobalScore_ch1(tp)
	    iC <- iC + 1
	}
}


# Plot histogram,
# superimpose histogram of predicted Global Scores
colRand <- "#E6E0FF"
colObs <- "#00DDAA44"
hist(scoreRand, 
     breaks=20, 
     xlim = c(min(min(scoreRand), min(scores$GLOBAL_SCORE)) * 1.2,
              max(max(scoreRand), max(scores$GLOBAL_SCORE)) * 1.2),
     freq = FALSE,
     col = colRand,
     main = "Random and observed Global Scores",
     cex.main = 0.8,
     xlab = "Global Score"
    )

if (length(scores$GLOBAL_SCORE) == 1) {
	abline(v = scores$GLOBAL_SCORE[rep(1,4)], col = colObs)
} else {
	par(new = TRUE)
	hVals <- hist(scores$GLOBAL_SCORE, plot = FALSE)
	hist(scores$GLOBAL_SCORE, 
	     xlim = c(min(min(scoreRand), min(scores$GLOBAL_SCORE)) * 1.2,
	              max(max(scoreRand), max(scores$GLOBAL_SCORE)) * 1.2),
	     ylim = c(0, 5 * max(hVals$density)),
	     breaks = round(N_RUNS)/2,
	     freq = FALSE,
	     col= colObs,
	     border= "#666677",
	     axes = FALSE,
	     main = "",
	     sub = (sprintf("mean Global Scores for predictions: %1.3f ± %1.3f",
	                    mean(scores$GLOBAL_SCORE), sd(scores$GLOBAL_SCORE))),
	     cex.sub = 0.8,
	     xlab = "",
	     ylab = ""
	    )
}

legend("topright",
       cex = 0.7,
       legend = c(sprintf("%d,%03d random",
                          floor(length(scoreRand)/1000),
                          length(scoreRand) %% 1000),
                  sprintf("%d predicted", N_RUNS)),
       fill = c(colRand, colObs),
       bty = "n",
       title = "")


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
     main = sprintf("%s crossvalidation with %d%% holdouts, %d runs\n(new correlations)",
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



# [END]
