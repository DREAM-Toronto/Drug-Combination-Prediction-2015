# analyzeSynergyDistributions.R
# simple script to evaluate ratio of synergetic
# experiments for a given combination. Some combinations
# are always synergetic, some never,  etc...

# Boris Jan. 28 


ABCSLOADER <- "loadABCSFromMaster.0.1.R" 
source(ABCSLOADER)

ABCS <- loadABCSmaster()
ABCS <- data.frame("ID" = character(nrow(ABCS)),
                          ABCS,
                          stringsAsFactors = FALSE)

ABCS$ID <- paste(ABCS$A, ABCS$B, sep=".")                  
ABCS <- ABCS[order(ABCS$ID),]

combi <- unique(ABCS$ID)

NC <- length(combi)
synDist <- data.frame("ID" = combi,        # combi ID
                      "S" = numeric(NC),   # number of synergetic experiments
                      "N" = numeric(NC),   # number of experiments
                      "Q" = numeric(NC),   # ratio S/N
                      stringsAsFactors = FALSE)

synThrsh <- 20   # threshold of synergy value to count an ABC combination as "synergetic"

# count number of "synergetic" experiments vs. all experiments
for (i in 1:length(combi)) {  # print combi-ID and synergetic vs. all
	syns <- ABCS$S[ABCS$ID == combi[i]]
	synDist$S[i] <- sum(as.numeric(syns > synThrsh)) 
	synDist$N[i] <- length(syns)
	synDist$Q[i] <- synDist$S[i] / synDist$N[i]
}

# order by ratio of S/N
synDist <- synDist[order(synDist$Q, decreasing = TRUE), ]


for (i in 1:nrow(synDist)) {  # print combi-ID and synergetic vs. all

	cat(sprintf("%40s: %d/%d\n",
	            synDist$ID[i],
	            synDist$S[i],
	            synDist$N[i]))
}


# [END]