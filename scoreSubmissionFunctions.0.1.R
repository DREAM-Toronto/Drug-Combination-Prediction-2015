# scoreSubmissionFunctions.0.1.R
#
# Purpose: Provides functions to score submissions for
#              cross-validation runs.
#          This file can be loaded as a SCORER asset in
#              the DREAM_Combi_main.R workflow and is
#              therefore expected to provide a
#              score() function that wraps/adapts the DREAM
#              challenge organizer's scoring scripts.
#
#          This file should contain only assets of functions
#             and constants. Sourcing this file must not have
#             any side-effects. Functions should not have side-
#             effects either.
#
# Version: 0.1           
#
# Date:    Jan 25 2016
# Author:  Boris and DREAM team UofT
#
# V 0.1    First code
#
# TODO:    
#
# ==========================================================

score <- function(type,
                  master,
                  pred,
                  conf,
                  topX = 100
                  ) {

	if (type == "1A" | type == "1B") {

	    # Default topX = 100 correlates all combination scores;
	    # Lower values select according to confidence rank.
		    
	    data <- SC$getObs(master, pred)
	    sc <- list()
        sc$CombiScores <- SC$getDrugCombiScore_ch1(data, conf, topX)
        sc$GlobalScore <- SC$getGlobalScore(data)

        return(sc)
        
	} else {
		stop("PANIC: type \"", type, "\" not yet supported.", sep="")
	}
}

# == INTERNAL FUNCTIONS ====================================

SC <- list()

SC$getObs <- function(master,
                      pred
                      ) {

    master <- data.frame(master,
                         "COMBINATION_ID" = sprintf("%s.%s", master$A, master$B),
                         stringsAsFactors=FALSE)
                         
	N <- nrow(pred)
	data <- data.frame("CELLS" = pred$CELL_LINE,
	                   "ID" = pred$COMBINATION_ID,
	                   "TRUE_SYN" = rep(as.numeric(NA), N),
	                   "PRED_SYN" = pred$PREDICTION,
	                   stringsAsFactors=FALSE)

    for (i in 1:nrow(data)) {
    	    syn  <- master[master$C == data$CELLS[i] &
    	                   master$COMBINATION_ID == data$ID[i],
    	                   "S"]
    	                   
    	    if (length(syn) == 1) {
    	    	    data$TRUE_SYN[i] <- syn
    	    } else if (length(syn) > 1) {
    	    	    data$TRUE_SYN[i] <- mean(syn)  # more than one observation
    	    	                                   # of this ABC combination
    	    	} else {
    	    		data$TRUE_SYN[i] <- NA
    	    	}          
    }

    if (any(is.na(data$TRUE_SYN))) {
    	    stop("PANIC: remaining NAs in TRUE_SYN column.")
    }
                                	
    return(data)                                	
}


SC$getDrugCombiScore_ch1 <- function(data, conf, topX) {
    # Calculate the drug combinations score of Subchallenge 1.
	# Modified from original code by Michael P Menden,
	#     Julio Saez-Rodriguez.
    
    # Also return MSE
    
    R <- matrix(numeric( 3 * nrow(conf)), ncol=3)
    colnames(R) <- c("R", "MSE", "CONFIDENCE")
    rownames(R) <- conf[, "COMBINATION_ID"]
    R[ , "CONFIDENCE"] <- conf[ , "CONFIDENCE"] 

    for (combID in rownames(R)) {
  	    true <- data[data[ , "ID"] == combID, "TRUE_SYN"]
  	    pred <- data[data[ , "ID"] == combID, "PRED_SYN"]
        R[combID, "R"] <- cor(true, pred)
# BS> This calculates correlations over _sets_ of 
# BS> combinations. This means the prediction is discarded when
# BS> a combination was used only once. cor is NA in that case.
# BS> When a combination is present twice, the correlation is
# BS> a perfect 1.0 or -1.0.
# BS> The leaderboard file contains 39/167 = 23.3% unique
# BS> combinations that are observed only twice! 
# BS> One should use a metric that gives robust results 
# BS> for small sample sizes.
        R[combID, "MSE"] <- mean((true - pred)^2)

     }

  #Make NA's in R <- 0
  R[is.na(R[ , "R"]), "R"] <- 0
# BS> This seems incorrect: one should remove NAs, not set them to 0
# BS> if one wants to calculate means. This will skew the mean unless
# BS> the mean IS zero.
  
  idx <- order(R[ , "CONFIDENCE"], decreasing = TRUE)
  idx <- idx[1:round(topX * (nrow(R) / 100))] # truncate to topX most
                                              # confident
  N <- sum(!is.na(R[idx, "R"]))
  return(c(meanR =   mean(R[idx, "R"],   na.rm=TRUE),
           steR  =     sd(R[idx, "R"],   na.rm=TRUE) / N,
           meanMSE = mean(R[idx, "MSE"], na.rm=TRUE),
           steMSE  =   sd(R[idx, "MSE"], na.rm=TRUE) / N,
           N    = N))
}


SC$getGlobalScore_ch1 <- function(data) {
    # Calculate the drug combinations GlobalScore of Subchallenge 1.
    #
	# Modified from original code by Michael P Menden,
	#     Julio Saez-Rodriguez.

    # Calculate Global Scores for subchallenge 1
    x <- data[, "TRUE_SYN"]
    y <- data[, "PRED_SYN"]
  
    agg <- aggregate(TRUE_SYN ~ CELLS, data, median)
    z0 <- agg$TRUE_SYN[match(data$CELLS, agg$CELLS)]
  
    agg <- aggregate(TRUE_SYN ~ ID, data, median)
    z1 <- agg$TRUE_SYN[match(data$ID, agg$ID)]
     
    # partial out the median of synergy across cell lines and combinations
    numerator <- SC$parCor_ch1(x,y,z1) - SC$parCor_ch1(x,z0,z1) * SC$parCor_ch1(z0,y,z1)
    denumerator <- sqrt(1-SC$parCor_ch1(x,z0,z1)^2) * sqrt(1-SC$parCor_ch1(z0,y,z1)^2)
  
    return(numerator/denumerator)
}


SC$parCor_ch1 <- function(u,v,w) {
    # partial correlations
    # originally by Michael P Menden, Julio Saez-Rodriguez
    numerator <- cor(u,v) - cor(u,w) * cor(w,v)
    denumerator <- sqrt(1-cor(u,w)^2) * sqrt(1-cor(w,v)^2)
    return(numerator/denumerator)
}


# [END]