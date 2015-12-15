setwd(DREAMDIR)

if (! require(RMySQL, quietly=TRUE)) {
  install.packages("RMySQL",repos="http://cran.us.r-project.org")
  library(RMySQL)
}

# == DATABASE LOGIN INFO ===========================
#
user <- 'zack'
pswd <- 'marchtenth0310'
dbname <- 'dream'
host <- 'localhost'

# == CONSTANTS =====================================
#
PRED_SAMPLE_DIR <- paste(DREAMDIR,"../Challenge Data/Drug Synergy Data/sample_submission",sep="")
CH1_LEAD <- "ch1_leaderboard-prediction.csv"
CH1_FINAL <- "ch1_final-prediction.csv"
CH2_LEAD <- "ch2_leaderboard-synergy_matrix.csv"
CH2_FINAL <- "ch2_final-synergy_matrix.csv"

# == HELPER FUNCTIONS ==============================
#

prepareMatrix <- function(triplets, features){
  preparedMatrix <- c()
  for(trip in triplets){
    cell <- trip[1]
    drug1 <- trip[2]
    drug2 <- trip[3]
    
    drug1Values <- drugValue(drug1,features.drug)
    drug2Values <- drugValue(drug2,features.drug)
    cellValues <- cellValue(cell,features.cell)
    drug1drug2Values <- drugDrugValue(drug1,drug2,features.drugDrug)
    drug1CellValues <- drugCellValue(drug1,cell,features.drugCell)
    drug2CellValues <- drugCellValue(drug2,cell,features.drugCell)
    
    tripValueList <- c(drug1Values,drug2Values,cellValues,drug1drug2Values,drug1CellValues,drug2CellValues)
    preparedMatrix <- c(preparedMatrix,tripValueList)
  }
  return(preparedMatrix)
}

drugDrugValue <- function(drug1,drug2,featureList){
  mydb = dbConnect(MySQL(), user=user, password=pswd, dbname=dbname, host=host)
  
  querySelect <- "select value from cell"
  queryJoin <- "join cell_name on cell_name.id=cell.cell_id join cell_feature on cell_feature.id=cell.feature_id"
  queryWhereCellName <- paste("cell_name.name= '",cell,"'",sep = "")
  
  # where cell_feature.feature=feature1 or cell_feature.feature=feature2 or .......
  queryWhereCellFeature <- paste("cell_feature.feature= '",featureList[1],"'",sep = "")
  if(length(featureList)>=2){
    for(feature in featureList[2:length(featureList)]){
      singleFeature <- paste("cell_feature.feature= '",feature,"'",sep = "")
      queryWhereCellFeature <- paste(queryWhereCellFeature,singleFeature,sep = " or ")
    }
  }
  
  queryWhere <- paste("where",queryWhereCellName,"and (",queryWhereCellFeature,")",sep = " ")
  # queryWhere <- "where cell_name.name='" + cell + "' and cell_feature.feature='" + 
  query <- paste(querySelect,queryJoin,queryWhere,sep = " ")
  
  rs = dbSendQuery(mydb,query)
  data = fetch(rs,n=-1)
  
  return(data)
}

cellValue <- function(cell,featureList){
  mydb = dbConnect(MySQL(), user=user, password=pswd, dbname=dbname, host=host)
  
  querySelect <- "select value from cell"
  queryJoin <- "join cell_name on cell_name.id=cell.cell_id join cell_feature on cell_feature.id=cell.feature_id"
  queryWhereCellName <- paste("cell_name.name= '",cell,"'",sep = "")
  
  # where cell_feature.feature=feature1 or cell_feature.feature=feature2 or .......
  queryWhereCellFeature <- paste("cell_feature.feature= '",featureList[1],"'",sep = "")
  if(length(featureList)>=2){
    for(feature in featureList[2:length(featureList)]){
      singleFeature <- paste("cell_feature.feature= '",feature,"'",sep = "")
      queryWhereCellFeature <- paste(queryWhereCellFeature,singleFeature,sep = " or ")
    }
  }
  
  queryWhere <- paste("where",queryWhereCellName,"and (",queryWhereCellFeature,")",sep = " ")
  # queryWhere <- "where cell_name.name='" + cell + "' and cell_feature.feature='" + 
  query <- paste(querySelect,queryJoin,queryWhere,sep = " ")
  
  rs = dbSendQuery(mydb,query)
  data = fetch(rs,n=-1)
  
  return(data)
}

drugValue <- function(drug,featureList){
  mydb = dbConnect(MySQL(), user=user, password=pswd, dbname=dbname, host=host)
  
  querySelect <- "select value from drug"
  queryJoin <- "join drug_name on drug_name.id=drug.drug_id join drug_feature on drug_feature.id=drug.feature_id"
  queryWhereDrugName <- paste("drug_name.name= '",drug,"'",sep = "")
  
  # where drug_feature.feature=feature1 or drug_feature.feature=feature2 or .......
  queryWhereDrugFeature <- paste("drug_feature.feature= '",featureList[1],"'",sep = "")
  if(length(featureList)>=2){
    for(feature in featureList[2:length(featureList)]){
      singleFeature <- paste("drug_feature.feature= '",feature,"'",sep = "")
      queryWhereDrugFeature <- paste(queryWhereDrugFeature,singleFeature,sep = " or ")
    }
  }
  
  queryWhere <- paste("where",queryWhereDrugName,"and (",queryWhereDrugFeature,")",sep = " ")
  # queryWhere <- "where drug_name.name='" + drug + "' and drug_feature.feature='" + 
  query <- paste(querySelect,queryJoin,queryWhere,sep = " ")
  
  rs = dbSendQuery(mydb,query)
  data = fetch(rs,n=-1)
  
  return(data)
}


getMatrix <- function(predictType,feature){
  chNum <- substr(predictType,1,1)
  chType <- substr(predictType,2,2)
  if(chNum=="1"){
    predFile <- (chType==l)? CH1_LEAD: CH1_FINAL
  } else{
    predFile <- (chType==l)? CH2_LEAD: CH2_FINAL
  }
}





findTriplet <- function(predType){
  predFile <- CH1_LEAD
  if(predType=='ch1f'){
    predFile <- CH1_TEST
  } else if(predType=='ch2l'){
    predFile <- CH2_LEAD
  } else if(predType=='ch2f'){
    predFile <- CH2_TEST
  } 
  predRawTemplate <- read.csv(paste(PRED_SAMPLE_DIR, predFile, sep=""),head=FALSE,stringsAsFactors=FALSE)
  predTripList <- predRawTemplate[2:nrow(predRawTemplate),1:3]
  
}