setwd(DREAMDIR)

if (! require(RMySQL, quietly=TRUE)) {
  install.packages("RMySQL",repos="http://cran.us.r-project.org")
  library(RMySQL)
}

# select drug_drug.value from drug_drug, drug_name as dnA, drug_name as dnB, drug_drug_feature where dnA.name='drugA' and dnB.name='drugC' and drug_drug_feature.feature='ddf2' and dnA.id = drug_drug.drugA_id and dnB.id=drug_drug.drugB_id and drug_drug.feature_id=drug_drug_feature.id;

# == DATABASE LOGIN INFO ===========================
#
user <- 'zack'
pswd <- 'marchtenth0310'
dbname <- 'dream'
host <- 'localhost'

# == CONSTANTS =====================================
#
PRED_SAMPLE_DIR <- paste(DREAMDIR,"sample_submission/",sep="")
CH1_LEAD <- "ch1_leaderboard-prediction.csv"
CH1_FINAL <- "ch1_final-prediction.csv"
CH2_LEAD <- "ch2_leaderboard-synergy_matrix.csv"
CH2_FINAL <- "ch2_final-synergy_matrix.csv"

# == HELPER FUNCTIONS ==============================
#

getMatrix <- function(chType, featureList){
  # chType (challenge type) can be '1f', '1l', '2f', '2l'
  if(substr(chType,1,1)=='1'){
    if(substr(chType,2,2)=='l'){
      chFile <- CH1_LEAD
    }else{
      chFile <- CH1_FINAL
    }
    ch1Sample <- read.csv(paste(PRED_SAMPLE_DIR, chFile, sep=""),head=FALSE,stringsAsFactors=FALSE)
    featureListLength <- length(featureList[[1]])+length(featureList[[2]])+length(featureList[[3]])+length(featureList[[4]])
    ch1Matrix <- matrix(nrow = (nrow(ch1Sample)*4),ncol = featureListLength)
    i <- 0
    for(triplet in ch1Sample[2:nrow(ch1Sample),]){
      i <- i+1
      cell <- triplet[1]
      drugPair <- strsplit(triplet[2],"[.]")
      drugA <- drugPair[1]
      drugB <- drugPair[2]
      tripletMatrix <- getSingleTripletMatrix(list(drugA,drugB,cell),featureList)
      ch1Matrix[i*4-3:i*4,] <- tripletMatrix
    }
    return(ch1Matrix)
    
    
    
  }else{
    if(substr(chType,2,2)=='l'){
      chFile <- CH2_LEAD
    }else{
      chFile <- CH2_FINAL
    }
    ch2Sample <- read.csv(paste(PRED_SAMPLE_DIR, chFile, sep=""),head=FALSE,stringsAsFactors=FALSE)
  }
}

getFeatureQuery <- function(features){
  featureQuery <- ""
  for(feature in features){
    if(featureQuery==""){
      featureQuery <- paste("and (f_table.feature='",feature,"'",sep = "")
    }else{
      single_feature <- paste("f_table.feature='",feature,"'",sep = "")
      featureQuery <- paste(featureQuery,single_feature,sep = " or ")
    }
  }
  featureQuery <- paste(featureQuery,");",sep = "")
  return(featureQuery)
}

getSingleTripletMatrix <- function(triplet, featureList){
  mydb = dbConnect(MySQL(), user=user, password=pswd, dbname=dbname, host=host)
  
  drugFeatures <- featureList[[1]]
  cellFeatures <- featureList[[2]]
  drugDrugFeatures <- featureList[[3]]
  drugCellFeatures <- featureList[[4]]
  # drugDrugCellFeatures <- featureList[[5]]
  
  preparedMatrix <- matrix(nrow = 4, ncol = (length(drugFeatures)+length(cellFeatures)+length(drugDrugFeatures)+length(drugCellFeatures)))
  
  drugA <- triplet[[1]]
  drugB <- triplet[[2]]
  cell <- triplet[[3]]
  
  drugAQuery1 <- paste("select drug.value from drug, drug_name as dnA, drug_feature as f_table where dnA.id = drug.drug_id and drug.feature_id=f_table.id and dnA.name='",drugA,"'",sep = "")
  drugAQuery2 <- getFeatureQuery(drugFeatures)
  drugAQuery <- paste(drugAQuery1,drugAQuery2,sep = " ")
  rs = dbSendQuery(mydb,drugAQuery)
  drugAData = fetch(rs,n=-1)
  drugAData <- t(as.matrix(drugAData))
  
  drugBQuery1 <- paste("select drug.value from drug, drug_name as dnA, drug_feature as f_table where dnA.id = drug.drug_id and drug.feature_id=f_table.id and dnA.name='",drugB,"'",sep = "")
  drugBQuery2 <- getFeatureQuery(drugFeatures)
  drugBQuery <- paste(drugBQuery1,drugBQuery2,sep = " ")
  rs = dbSendQuery(mydb,drugBQuery)
  drugBData = fetch(rs,n=-1)
  drugBData <- t(as.matrix(drugBData))
  
  cellQuery1 <- paste("select cell.value from cell, cell_name as cn, cell_feature as f_table where cn.id = cell.cell_id and cell.feature_id=f_table.id and cn.name='",cell,"'",sep = "")
  cellQuery2 <- getFeatureQuery(cellFeatures)
  cellQuery <- paste(cellQuery1,cellQuery2,sep = " ")
  rs = dbSendQuery(mydb,cellQuery)
  cellData = fetch(rs,n=-1)
  cellData <- t(as.matrix(cellData))
  
  drugDrugQuery1 <- paste("select drug_drug.value from drug_drug, drug_name as dnA, drug_name as dnB, drug_drug_feature as f_table where dnA.id = drug_drug.drugA_id and dnB.id = drug_drug.drugB_id and drug_drug.feature_id=f_table.id and dnA.name='",drugA,"' and dnB.name='",drugB,"'",sep = "")
  drugDrugQuery2 <- getFeatureQuery(drugDrugFeatures)
  drugDrugQuery <- paste(drugDrugQuery1,drugDrugQuery2,sep = " ")
  rs = dbSendQuery(mydb,drugDrugQuery)
  drugDrugData = fetch(rs,n=-1)
  drugDrugData <- t(as.matrix(drugDrugData))
  
  drugACellQuery1 <- paste("select drug_cell.value from drug_cell, drug_name as dnA, cell_name as cn, drug_cell_feature as f_table where dnA.id = drug_cell.drug_id and cn.id = drug_cell.cell_id and drug_cell.feature_id=f_table.id and dnA.name='",drugA,"' and cn.name='",cell,"'",sep = "")
  drugACellQuery2 <- getFeatureQuery(drugCellFeatures)
  drugACellQuery <- paste(drugACellQuery1,drugACellQuery2,sep = " ")
  rs = dbSendQuery(mydb,drugACellQuery)
  drugACellData = fetch(rs,n=-1)
  drugACellData <- t(as.matrix(drugACellData))
  
  drugBCellQuery1 <- paste("select drug_cell.value from drug_cell, drug_name as dnA, cell_name as cn, drug_cell_feature as f_table where dnA.id = drug_cell.drug_id and cn.id = drug_cell.cell_id and drug_cell.feature_id=f_table.id and dnA.name='",drugB,"' and cn.name='",cell,"'",sep = "")
  drugBCellQuery2 <- getFeatureQuery(drugCellFeatures)
  drugBCellQuery <- paste(drugBCellQuery1,drugBCellQuery2,sep = " ")
  rs = dbSendQuery(mydb,drugBCellQuery)
  drugBCellData = fetch(rs,n=-1)
  drugBCellData <- t(as.matrix(drugBCellData))

    # drugAData, drugBData, cellData, drugDrugData, drugACellData, drugBCellData, drugDrugCellData
    abab <- cbind(drugAData, drugBData, cellData, drugDrugData, drugACellData, drugBCellData)
    abba <- cbind(drugAData, drugBData, cellData, drugDrugData, drugBCellData, drugACellData)
    baab <- cbind(drugBData, drugAData, cellData, drugDrugData, drugACellData, drugBCellData)
    baba <- cbind(drugBData, drugAData, cellData, drugDrugData, drugBCellData, drugACellData)
    
#   # drugAData, drugBData, cellData, drugDrugData, drugACellData, drugBCellData, drugDrugCellData
#   preparedMatrix[1,] <- cbind(drugAData, drugBData, cellData, drugDrugData, drugACellData, drugBCellData)
#   preparedMatrix[2,] <- cbind(drugAData, drugBData, cellData, drugDrugData, drugBCellData, drugACellData)
#   preparedMatrix[3,] <- cbind(drugBData, drugAData, cellData, drugDrugData, drugACellData, drugBCellData)
#   preparedMatrix[4,] <- cbind(drugBData, drugAData, cellData, drugDrugData, drugBCellData, drugACellData)
  
  currTripValues <- rbind(abab,abba,baab,baba)
  
  return(currTripValues)
}