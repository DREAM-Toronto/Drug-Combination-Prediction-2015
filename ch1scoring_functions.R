# ch1scoring_functions.R
# These are the AZ Challenge organizers' scoring functions
# modified to work as local functions in xValidation
# runs.
#
# Version: 1.0           
#
# Date:    Nov 9 2015
# Author:  Michael P Menden, Julio Saez-Rodriguez;
#          modifications by Boris

# V 1.0    Adapted to work as functions in our workflow
# ==========================================================

# == CONSTANTS =============================================
#

# == FUNCTIONS =============================================
#


getDrugCombiScore_ch1 <- function(data, conf, topX = 100) {
    # Calculate the drug combinations score of Subchallenge 1.
	# Modified from original code by Michael P Menden,
	#     Julio Saez-Rodriguez.
    # Default topX = 100 correlates all combination scores;
    # lower values correlate according to confidence in
    # column SYN_CONF
    
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


parCor <- function(u,v,w) {
    # partial correlations
    # originally by Michael P Menden, Julio Saez-Rodriguez
    numerator <- cor(u,v) - cor(u,w) * cor(w,v)
    denumerator <- sqrt(1-cor(u,w)^2) * sqrt(1-cor(w,v)^2)
    return(numerator/denumerator)
}



getGlobalScore_ch1 <- function(data) {
    # Calculate Global Scores for subchallenge 1
    x <- data[, "TRUE_SYN"]
    y <- data[, "PRED_SYN"]
  
    agg <- aggregate(TRUE_SYN ~ CELLS, data, median)
    z0 <- agg$TRUE_SYN[match(data$CELLS, agg$CELLS)]
  
    agg <- aggregate(TRUE_SYN ~ ID, data, median)
    z1 <- agg$TRUE_SYN[match(data$ID, agg$ID)]
     
    # partial out the median of synergy across cell lines and combinations
    numerator <- parCor(x,y,z1) - parCor(x,z0,z1) * parCor(z0,y,z1)
    denumerator <- sqrt(1-parCor(x,z0,z1)^2) * sqrt(1-parCor(z0,y,z1)^2)
  
    return(numerator/denumerator)
}






