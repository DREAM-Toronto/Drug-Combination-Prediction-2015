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

a <- list(b<-3,c<-4)
d <- 0
for(i in a){
  d<-d+i
  print(d)
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



# select drug_drug.value from drug_drug, drug_name as dnA, drug_name as dnB, drug_drug_feature where dnA.name='drugA' and dnB.name='drugC' and drug_drug_feature.feature='ddf2' and dnA.id = drug_drug.drugA_id and dnB.id=drug_drug.drugB_id and drug_drug.feature_id=drug_drug_feature.id;
# tripletList[[i]] -> one triplet and tripletList[[i]][1] -> cell, tripletList[[i]][2] -> drugA..
# featureList[[1]] -> drug feature list and featureList[[1]][[1]] -> first feature in drug feature list
prepareMatrix <- function(tripletList, featureList){
  drugFeatures <- featureList[[1]]
  cellFeatures <- featureList[[2]]
  drugDrugFeatures <- featureList[[3]]
  drugCellFeatures <- featureList[[4]]
  drugDrugCellFeatures <- featureList[[5]]
  
  index <- 0
  
  preparedMatrix <- matrix(nrow = nrow(tripletList),ncol = (length(drugFeatures)+length(cellFeatures)+length(drugDrugFeatures)+length(drugCellFeatures)))
  
  for(trip in triplets){
    index <- index + 1
    drugA <- trip[[2]]
    drubB <- trip[[3]]
    cell <- trip[[1]]
    
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
    
#     drugDrugCellQuery1 <- paste("select drug_drug_cell.value from drug_drug_cell, drug_name as dnA, drug_name as dnB, cell_name as cn, drug_drug_cell_feature as f_table where dnA.id = drug_drug_cell.drugA_id and dnB.id = drug_drug_cell.drugB_id and cn.id = drug_drug_cell.cell_id and drug_drug_cell.feature_id=f_table.id and dnA.name='",drugA,"' and dnB.name='",drugB,"' and cn.name='",cell,"'",sep = "")
#     drugDrugCellQuery2 <- getFeatureQuery(features)
#     drugDrugCellQuery <- paste(drugDrugCellQuery1,drugDrugCellQuery2,sep = " ")
#     rs = dbSendQuery(mydb,drugDrugCellQuery)
#     drugDrugCellData = fetch(rs,n=-1)
#     drugDrugCellData <- t(as.matrix(drugDrugCellData))
#     
#     # drugAData, drugBData, cellData, drugDrugData, drugACellData, drugBCellData, drugDrugCellData
#     abab <- cbind(drugAData, drugBData, cellData, drugDrugData, drugACellData, drugBCellData, drugDrugCellData)
#     abba <- cbind(drugAData, drugBData, cellData, drugDrugData, drugBCellData, drugACellData, drugDrugCellData)
#     baab <- cbind(drugBData, drugAData, cellData, drugDrugData, drugACellData, drugBCellData, drugDrugCellData)
#     baba <- cbind(drugBData, drugAData, cellData, drugDrugData, drugBCellData, drugACellData, drugDrugCellData)
    
    currTripValues <- rbind(abab,abba,baab,baba)
    
    preparedMatrix[i*4-3:i*4,] <- currTripValues
  }
}



prepareMatrixOLD <- function(triplets, features){
  mydb = dbConnect(MySQL(), user=user, password=pswd, dbname=dbname, host=host)
  # preparedMatrix <- c()
  for(trip in triplets){
    drugA <- trip[2]
    drubB <- trip[3]
    cell <- trip[1]
    
    drugQuery1 <- paste("select drug.value from drug, drug_name as dnA, drug_feature as f_table where dnA.id = drug.drug_id and drug.feature_id=f_table.id and dnA.name='",drugA,"'",sep = "")
    drugQuery2 <- getFeatureQuery(features)
    drugQuery <- paste(drugQuery1,drugQuery2,sep = " ")
    rs = dbSendQuery(mydb,drugQuery)
    drugData = fetch(rs,n=-1)
    
    cellQuery1 <- paste("select cell.value from cell, cell_name as cn, cell_feature as f_table where cn.id = cell.cell_id and cell.feature_id=f_table.id and cn.name='",cell,"'",sep = "")
    cellQuery2 <- getFeatureQuery(features)
    cellQuery <- paste(cellQuery1,cellQuery2,sep = " ")
    rs = dbSendQuery(mydb,cellQuery)
    cellData = fetch(rs,n=-1)
    
    drugDrugQuery1 <- paste("select drug_drug.value from drug_drug, drug_name as dnA, drug_name as dnB, drug_drug_feature as f_table where dnA.id = drug_drug.drugA_id and dnB.id = drug_drug.drugB_id and drug_drug.feature_id=f_table.id and dnA.name='",drugA,"' and dnB.name='",drugB,"'",sep = "")
    drugDrugQuery2 <- getFeatureQuery(features)
    drugDrugQuery <- paste(drugDrugQuery1,drugDrugQuery2,sep = " ")
    rs = dbSendQuery(mydb,drugDrugQuery)
    drugDrugData = fetch(rs,n=-1)
    
    drugCellQuery1 <- paste("select drug_cell.value from drug_cell, drug_name as dnA, cell_name as cn, drug_cell_feature as f_table where dnA.id = drug_cell.drug_id and cn.id = drug_cell.cell_id and drug_cell.feature_id=f_table.id and dnA.name='",drugA,"' and cn.name='",cell,"'",sep = "")
    drugCellQuery2 <- getFeatureQuery(features)
    drugCellQuery <- paste(drugCellQuery1,drugCellQuery2,sep = " ")
    rs = dbSendQuery(mydb,drugCellQuery)
    drugCellData = fetch(rs,n=-1)
    
    drugDrugCellQuery1 <- paste("select drug_drug_cell.value from drug_drug_cell, drug_name as dnA, drug_name as dnB, cell_name as cn, drug_drug_cell_feature as f_table where dnA.id = drug_drug_cell.drugA_id and dnB.id = drug_drug_cell.drugB_id and cn.id = drug_drug_cell.cell_id and drug_drug_cell.feature_id=f_table.id and dnA.name='",drugA,"' and dnB.name='",drugB,"' and cn.name='",cell,"'",sep = "")
    drugDrugCellQuery2 <- getFeatureQuery(features)
    drugDrugCellQuery <- paste(drugDrugCellQuery1,drugDrugCellQuery2,sep = " ")
    rs = dbSendQuery(mydb,drugDrugCellQuery)
    drugDrugCellData = fetch(rs,n=-1)
    
    drugQuery <- "select drug_drug.value from drug_drug, drug_name as dnA, drug_name as dnB, drug_drug_feature where dnA.name='drugA' and dnB.name='drugC' and drug_drug_feature.feature='ddf2' and dnA.id = drug_drug.drugA_id and dnB.id=drug_drug.drugB_id and drug_drug.feature_id=drug_drug_feature.id;"
    
#     drug1 <- trip[2]
#     drug2 <- trip[3]
#     cell <- trip[1]
#     
#     drug1Values <- drugValue(drug1,features.drug)
#     drug2Values <- drugValue(drug2,features.drug)
#     cellValues <- cellValue(cell,features.cell)
#     drug1drug2Values <- drugDrugValue(drug1,drug2,features.drugDrug)
#     drug1CellValues <- drugCellValue(drug1,cell,features.drugCell)
#     drug2CellValues <- drugCellValue(drug2,cell,features.drugCell)
#     
#     tripValueList <- c(drug1Values,drug2Values,cellValues,drug1drug2Values,drug1CellValues,drug2CellValues)
#     preparedMatrix <- c(preparedMatrix,tripValueList)
  }
  return(preparedMatrix)
}


# Need to complete drugCellValue
drugCellValue <- function(drug,cell,featureList){
  mydb = dbConnect(MySQL(), user=user, password=pswd, dbname=dbname, host=host)
  
  querySelect <- "select value from drug_cell"
  queryJoin1 <- "join drug_name on drug_cell.drugA_id=drugA_name.id";
  queryJoin2 <- "join drug_name as drugB_name on drug_drug.drugB_id=drugB_name.id";
  queryJoin3 <- "join drug_feature on drug_drug.feature_id=drug_feature.id";
  queryJoin <- paste(queryJoin1,queryJoin2,queryJoin3,sep = " ");
  
  queryWhereNames <- paste("where drugA_name.name='",drug1,"' and drugB_name.name='",drug2,"' and (",sep = "")
  queryWhereFeature <- ""
  for(feature in featureList){
    if(queryWhereFeature==""){
      queryWhereFeature <- paste("drug_drug_feature.feature='",feature,"'",sep = "")
    }else{
      singleFeature <- paste("drug_drug_feature.feature='",feature,"'",sep = "")
      queryWhereFeature <- paste(queryWhereFeature,singleFeature,sep = " or ")  
    }
  }
  queryWhere <- paste(queryWhereNames,queryWhereFeature,")",sep = "")
  
  query <- paste(querySelect,queryJoin,queryWhere)
  
  
  rs = dbSendQuery(mydb,query)
  data = fetch(rs,n=-1)
  
  return(data)
}

drugDrugValue <- function(drug1,drug2,featureList){
  mydb = dbConnect(MySQL(), user=user, password=pswd, dbname=dbname, host=host)
  
  querySelect <- "select value from drug_drug"
  queryJoin1 <- "join drug_name as drugA_name on drug_drug.drugA_id=drugA_name.id";
  queryJoin2 <- "join drug_name as drugB_name on drug_drug.drugB_id=drugB_name.id";
  queryJoin3 <- "join drug_feature on drug_drug.feature_id=drug_feature.id";
  queryJoin <- paste(queryJoin1,queryJoin2,queryJoin3,sep = " ");
  
  queryWhereNames <- paste("where drugA_name.name='",drug1,"' and drugB_name.name='",drug2,"' and (",sep = "")
  queryWhereFeature <- ""
  for(feature in featureList){
    if(queryWhereFeature==""){
      queryWhereFeature <- paste("drug_drug_feature.feature='",feature,"'",sep = "")
    }else{
      singleFeature <- paste("drug_drug_feature.feature='",feature,"'",sep = "")
      queryWhereFeature <- paste(queryWhereFeature,singleFeature,sep = " or ")  
    }
  }
  queryWhere <- paste(queryWhereNames,queryWhereFeature,")",sep = "")
  
  query <- paste(querySelect,queryJoin,queryWhere)
  # return(query)
  
#   queryWhereCellName <- paste("cell_name.name= '",cell,"'",sep = "")
#   
#   # where cell_feature.feature=feature1 or cell_feature.feature=feature2 or .......
#   queryWhereCellFeature <- paste("cell_feature.feature= '",featureList[1],"'",sep = "")
#   if(length(featureList)>=2){
#     for(feature in featureList[2:length(featureList)]){
#       singleFeature <- paste("cell_feature.feature= '",feature,"'",sep = "")
#       queryWhereCellFeature <- paste(queryWhereCellFeature,singleFeature,sep = " or ")
#     }
#   }
#   
#   queryWhere <- paste("where",queryWhereCellName,"and (",queryWhereCellFeature,")",sep = " ")
#   # queryWhere <- "where cell_name.name='" + cell + "' and cell_feature.feature='" + 
#   query <- paste(querySelect,queryJoin,queryWhere,sep = " ")
  
#   query1 <- "select value from drug_drug join drug_name as drugA_name on drug_drug.drugA_id=drugA_name.id"
#   query2 <- "join drug_name as drugB_name on drug_drug.drugB_id=drugB_name.id"
#   query3 <- "join drug_feature on drug_drug.feature_id=drug_feature.id where drugA_name.name='"
#   query4 <- drug1
#   query5 <- "' and drugB_name.name='"
#   query6 <- drug2
#   query7 <- "' and "
#   query7 <- "';"
#   query8 <- 
  # query <- "select value from drug_drug join drug_name as drugA_name on drug_drug.drugA_id=drugA_name.id join drug_name as drugB_name on drug_drug.drugB_id=drugB_name.id join drug_feature on drug_drug.feature_id=drug_feature.id where drugA_name.name='drugB' and drugB_name.name='drugA';"
  
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

findCh1Triplet <- function(fn){
  predRawTemplate <- read.csv(paste(PRED_SAMPLE_DIR, fn, sep=""),head=FALSE,stringsAsFactors=FALSE)
  predTripList <- predRawTemplate[2:nrow(predRawTemplate),1:2]
}



findTriplet <- function(predType){
  predFile <- ""
  if(substr(predType,1,3)=='ch1'){
    if(substr(predType,4,4)=='l'){
      predFile <- CH1_LEAD
    }else{
      predFile <- CH1_FINAL
    }
    return(findCh1Triplet(predFile))
  }else{
    if(substr(predType,4,4)=='l'){
      predFile <- CH2_LEAD
    }else{
      predFile <- CH2_FINAL
    }
    findCh2Triplet(predFile)
  }
}


findTripletOLD <- function(predType){
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