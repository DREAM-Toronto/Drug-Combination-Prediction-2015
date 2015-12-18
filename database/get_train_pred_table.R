setwd(DREAMDIR)

getCh1TraPreTriplets <- function(){
  trainPath <- paste(DREAMDIR,"../Challenge Data/Drug Synergy Data/ch1_train_combination_and_monoTherapy.csv",sep = "")
  ch1LeadPath <- paste(DREAMDIR,"sample_submission/ch1_leaderboard-prediction.csv",sep = "")
  trainData <- read.csv(trainPath,header = FALSE,stringsAsFactors=FALSE)
  ch1LeadData <- read.csv(ch1LeadPath,header = FALSE,stringsAsFactors=FALSE)
  
  trainTriplets <- as.matrix(trainData[2:nrow(trainData),1:3])
  trainTriplets[,c(1,2,3)] <- trainTriplets[,c(2,3,1)]
  # ch1LeadTriplets
  
  ch1Cell <- matrix(ch1LeadData[2:nrow(ch1LeadData),1],ncol = 1)
  ch1DrugPairList <- t(matrix(unlist(strsplit(ch1LeadData[2:nrow(ch1LeadData),2],"[.]")),nrow = 2))
  ch1LeadTriplets <- cbind(ch1DrugPairList,ch1Cell)
  
  ch1TraPreTriplets <- rbind(trainTriplets,ch1LeadTriplets)
  rowNames <- 1:nrow(ch1TraPreTriplets)
  ch1TraPreTriplets <- data.frame(ch1TraPreTriplets,row.names = rowNames)
  names(ch1TraPreTriplets) <- c("drugA","drugB","cell")
  return(ch1TraPreTriplets)
}
