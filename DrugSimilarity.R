MONO_FILE <- "../Challenge Data/Drug Synergy Data/Drug_info_release.csv"
DRUG_PROPERTIES <- "../Challenge Data/Drug Synergy Data/Drug_info_release.csv"
CELL_KILL_DIR <- "../Challenge Data/Drug Synergy Data/Raw Data/Raw_Data_csv/ch1_ch2_monoTherapy/"
CELL_PROPERTIES <- "../Challenge Data/Sanger Molecular Data/cell_info.csv"

# =============================================
# Should consider moving these functions into utilities
# Get the list all unique drug names 
getAllDrugs <- function() {
  raw <- read.csv(paste("../Challenge Data/Drug Synergy Data/", "Drug_info_release.csv", sep=""),
                  head=FALSE,
                  stringsAsFactors=FALSE)
  return(raw[2:dim(raw)[1],1])
}

# Get the list of all unique cell lines 
getAllCellLines <- function() {
  raw <- read.csv(paste("../Challenge Data/Sanger Molecular Data/", "cell_info.csv", sep=""),
                  head=FALSE,
                  stringsAsFactors=FALSE)
  return(raw[2:dim(raw)[1],1])
  
}

# drug is drug object containing s
getCurveParams <- function (drug, cellLine) {
  
  dose <- seq (0, 1, 0.01)
  
  # stores into a data frame, mono therapy effects given a specific drug and cell line
  dfMonoEffect <- getMonoEffect(drug, cellLine)
  
  # get drug response curves 
  drc <- list()
  # normalize the concentrations
  drc$conc <- c(0, 0.01, 0.03, 0.1, 0.3, 1)
  # take the average of the effects 
  drc$dat <- colMeans(dfMonoEffect[8:13])
  
  # curve fit 
  coefDRC <- nlsDRC(drc)
  
  return(coefDRC)  
  
}


# =============================================
getMonoEffect <- function(drug,cellLine) {
  CELL_KILL_DIR <- "../Challenge Data/Drug Synergy Data/Raw Data/Raw_Data_csv/ch1_ch2_monoTherapy/"
  regEx <- paste(".*(",drug,"[.]).*(",cellLine,")",sep="")
  filenameList <- list.files(CELL_KILL_DIR,pattern = regEx)
  
  filename <- c()
  conc <- c()
  eff <- c()
  
  for(fn in filenameList){
    raw <- read.csv(paste(CELL_KILL_DIR, fn, sep=""),
                    head=FALSE,
                    stringsAsFactors=FALSE)
    
    filename <- c(filename,fn)
    
    if(raw[9,2]==drug){
      # Agent1
      conc <- as.numeric(c(conc,raw[2:7,1]))
      eff <- as.numeric(c(eff,raw[2:7,2]))
    }else{
      # Agent2
      conc <- as.numeric(c(conc,raw[1,2:7]))
      eff <- as.numeric(c(eff,raw[2,2:7]))
    }
  }
  
  concentration <- matrix(conc,length(filenameList),6,TRUE)
  effectiveness <- matrix(eff,length(filenameList),6,TRUE)
  
  allMonoEffects <- data.frame(filename,concentration,effectiveness)
  
  return(allMonoEffects)
}

# == PARSE DRUG INFO CSV ========================================================
# Get chemical/physical properties of a drug

getDrugProperties <- function(drug){
  raw <- read.csv(MONO_FILE,head=FALSE,stringsAsFactors=FALSE)
  
  nDrugs <- (dim(raw))[1]-1
  drugList <- raw[,1]
  drugLoc <- match(drug,drugList)
  
  rawDrugProperties <- raw[drugLoc,2:8]
  
  drugProperties <- list()
  
  # Store member targets as a vector
  rawDrugTargetInString <- rawDrugProperties[1,1]
  rawDrugTargetInString <- as.character(rawDrugTargetInString)
  rawDrugTarget <- strsplit(rawDrugTargetInString,",")
  rawDrugTargetInVector <- c(rawDrugTarget)
  drugProperties$Target <- rawDrugTargetInVector
  
  # Store member HBA
  drugProperties$HBA <- as.numeric(rawDrugProperties[1,2])
  
  # Store member CLogP
  drugProperties$CLogP <- as.numeric(rawDrugProperties[1,3])
  
  # Store member HBD
  drugProperties$HBD <- as.numeric(rawDrugProperties[1,4])
  
  # Store member MW
  drugProperties$MW <- as.numeric(rawDrugProperties[1,8])
  
  return(drugProperties)
}


# == DISTANCE METRICS ========================================================
# Defining distance metrics for each feature of a drug

distTarget <- function(target1, target2) {
  # modified jaccard distance
  # taking into account of number of targets given
  diff <- 1 - length(intersect(target1, target2))
  return (diff^2/length(union(target1, target2)))
}

distHBA <- function(HBA1, HBA2) {
  abs(HBA1 - HBA2)
}

distHBD <- function(HBD1, HBD2) {
  abs(HBD1 - HBD2)
}

distCLogP <- function(CLogP1, CLogP2) {
  10^abs(CLogP1 - CLogP2)
}

distMW <- function(MW1, MW2) {
  abs(MW1 - MW2)  
}

# compute the distance between two drugs
distDrug <- function(drug1, drug2) {
  
  distance <- c(
    distTarget(drug1[[1]],drug2[[1]]),
    distHBA(drug1[[2]],drug2[[2]]),
    distCLogP(drug1[[3]],drug2[[3]]),
    distHBD(drug1[[4]],drug2[[4]]),
    distMW(drug1[[5]],drug2[[5]])
  )
  return(distance)
}

# == GETTING FEATURE MATRICES FROM MONO THERAPY DATA ============================
getDrugFeatureMatrix <- function(){
  monoFeatureMatrix <- getMonoFeatureMatrix()
  fillMonoFeatureWithCellAverage(monoFeatureMatrix)
}

getCellFeatureMatrix <- function(){
  monoFeatureMatrix <- getMonoFeatureMatrix()
  cfm <- fillMonoFeatureWithDrugAverage(monoFeatureMatrix)
  return (t(cfm))
}

getMonoFeatureMapFromFiles <- function(filePathVect){
  monoFeatureMap <- list()
  for(filePath in filePathVect){
    monoFeatureRaw <- read.csv(filePath,
                               head=FALSE,
                               stringsAsFactors=FALSE)
    nRow <- dim(monoFeatureRaw)[1]
    for(i in 2:nRow){
      cellLine <- monoFeatureRaw[i,1]
      
      # For drug1
      drug1 <- monoFeatureRaw[i,2]
      drug1CellKey <- paste(drug1,",",cellLine,sep="")
      
      maxConcentration1 <- as.numeric(monoFeatureRaw[i,4])
      
      if(is.null(monoFeatureMap[[drug1CellKey]])){
        ic50<-c()
        h<-c()
        einf<-c()
        
        drugCellValue<-list()
        drugCellValue$IC50 <- ic50
        drugCellValue$H <- h
        drugCellValue$Einf <- einf
        
        monoFeatureMap[[drug1CellKey]] <- drugCellValue
      }
      
      monoFeatureMap[[drug1CellKey]]$IC50 <- c(monoFeatureMap[[drug1CellKey]]$IC50,(as.numeric(monoFeatureRaw[i,6]))/maxConcentration1)
      monoFeatureMap[[drug1CellKey]]$H <- c(monoFeatureMap[[drug1CellKey]]$H,as.numeric(monoFeatureRaw[i,7])*maxConcentration1)
      monoFeatureMap[[drug1CellKey]]$Einf <- c(monoFeatureMap[[drug1CellKey]]$Einf,as.numeric(monoFeatureRaw[i,8]))
      
      # For drug2
      drug2 <- monoFeatureRaw[i,3]
      drug2CellKey <- paste(drug2,",",cellLine,sep="")
      
      maxConcentration2 = as.numeric(monoFeatureRaw[i,5])
      
      if(is.null(monoFeatureMap[[drug2CellKey]])){
        ic50<-c()
        h<-c()
        einf<-c()
        
        drugCellValue<-list()
        drugCellValue$IC50 <- ic50
        drugCellValue$H <- h
        drugCellValue$Einf <- einf
        
        monoFeatureMap[[drug2CellKey]] <- drugCellValue
      }
      
      monoFeatureMap[[drug2CellKey]]$IC50 <- c(monoFeatureMap[[drug2CellKey]]$IC50,(as.numeric(monoFeatureRaw[i,9]))/maxConcentration2)
      monoFeatureMap[[drug2CellKey]]$H <- c(monoFeatureMap[[drug2CellKey]]$H,as.numeric(monoFeatureRaw[i,10])*maxConcentration2)
      monoFeatureMap[[drug2CellKey]]$Einf <- c(monoFeatureMap[[drug2CellKey]]$Einf,as.numeric(monoFeatureRaw[i,11]))
    }
  }
  nUniquePair <- length(monoFeatureMap)
  for(i in 1:nUniquePair){
    drugCellValue <- monoFeatureMap[[i]]
    ic50 <- mean(drugCellValue$IC50)
    h <- mean(drugCellValue$H)
    einf <- mean(drugCellValue$Einf)
    monoFeatureMap[[i]]$IC50 <- ic50
    monoFeatureMap[[i]]$H <- h
    monoFeatureMap[[i]]$Einf <- einf
  }
  return(monoFeatureMap)
}

getMonoFeatureMatrix <- function(){
  pathVect <- c("../Challenge Data/Drug Synergy Data/ch1_train_combination_and_monoTherapy.csv",
                "../Challenge Data/Drug Synergy Data/ch1_leaderBoard_monoTherapy.csv",
                "../Challenge Data/Drug Synergy Data/ch1_test_monoTherapy.csv")
  monoFeatureMap <- getMonoFeatureMapFromFiles(pathVect)
  
  drugListFilePath <- "../Challenge Data/Drug Synergy Data/Drug_info_release.csv"
  raw <- read.csv(drugListFilePath,
                  head=FALSE,
                  stringsAsFactors=FALSE)
  drugVect <- raw[-1,1]
  
  cellListFilePath <- "../Challenge Data/Sanger Molecular Data/cell_info.csv"
  raw <- read.csv(cellListFilePath,
                  head=FALSE,
                  stringsAsFactors=FALSE)
  cellVect <- raw[-1,1]
  
  monoFeatureMatrix <- matrix(data=NA,nrow=length(drugVect),ncol=3*length(cellVect))
  
  for(i in 1:length(drugVect)){
    for(j in 1:length(cellVect)){
      drug <- drugVect[i]
      cellLine <- cellVect[j]
      drugCellKey <- paste(drug,",",cellLine,sep="")
      if(!is.null(monoFeatureMap[[drugCellKey]])){
        monoFeatureMatrix[i,j*3-2] <- monoFeatureMap[[drugCellKey]]$IC50
        monoFeatureMatrix[i,j*3-1] <- monoFeatureMap[[drugCellKey]]$H
        monoFeatureMatrix[i,j*3] <- monoFeatureMap[[drugCellKey]]$Einf
      }
    }
  }
  return(monoFeatureMatrix)
}


fillMonoFeatureWithCellAverage <- function(monoFeatureMatrix){
  nColTotal <- ncol(monoFeatureMatrix)
  nColForEachFeature <- nColTotal%/%3
  einfIndexSequence <- seq(3,nColTotal,3)
  hIndexSequence <- einfIndexSequence-1
  ic50IndexSequence <- einfIndexSequence-2
  
  ic50Average <- rowMeans(monoFeatureMatrix[,ic50IndexSequence],TRUE)
  hAverage <- rowMeans(monoFeatureMatrix[,hIndexSequence],TRUE)
  einfAverage <- rowMeans(monoFeatureMatrix[,einfIndexSequence],TRUE)
  cellAverage <- c(ic50Average,hAverage,einfAverage)
  
  missingDataLoc <- which(is.na(monoFeatureMatrix))
  missingDataReplacementIndexFromAverage <- missingDataLoc%%length(cellAverage)
  missingDataReplacementIndexFromAverage[missingDataReplacementIndexFromAverage==0] <- length(cellAverage)
  
  monoFeatureMatrix[missingDataLoc] <- cellAverage[missingDataReplacementIndexFromAverage]
  return(monoFeatureMatrix)
}

fillMonoFeatureWithDrugAverage <- function(monoFeatureMatrix){
  nRow <- nrow(monoFeatureMatrix)
  nCol <- ncol(monoFeatureMatrix)
  drugAverage <- colMeans(monoFeatureMatrix,TRUE)
  missingDataLoc <- which(is.na(monoFeatureMatrix))
  missingDataReplacementIndexFromAverage <- (missingDataLoc + nRow - 1) %/% nRow
  monoFeatureMatrix[missingDataLoc] <- drugAverage[missingDataReplacementIndexFromAverage]
  return(monoFeatureMatrix)
}

# == COMPUTE CORRELATION FROM MONO THERAPY DATA =======================================
# calculates correlation between two drugs (or two cell lines)
# return a vector of correlation by each cell line (by each drug)
# (the last value is the combined correlation using data from all cell lines (from all drugs))


corCellLine <- function (cell1, cell2) {
  
  corScores <- c()
  
  allCells <- getAllCellLines()
  
  dose <- seq (0, 1, 0.01)
  
  # cfm <- cellFeatureMatrix
  cfm <- getCellFeatureMatrix()
  
  start <- match(cell1, allCells)
  end <- match(cell1, allCells) + 2
  coefDRC1 <- cfm[start : end,]
  coefDRC2 <- cfm[start : end,]
  
  effect1 <- c()
  effect2 <- c()
  
  # for all cell lines
  for (i in seq(0,length(coefDRC1)/3 - 1)) {
    
    # interpolate data points from the curves
    newEffects1 <- fHill(dose, coefDRC1[1,i], coefDRC1[2,i], coefDRC1[3,i])
    newEffects2 <- fHill(dose, coefDRC2[1,i], coefDRC2[2,i], coefDRC2[3,i])
    effect1 <- c(effect1, newEffects1)
    effect2 <- c(effect2, newEffects2)
    
    corScores <- c(corScores, cor(newEffects1, newEffects2, "everything", "pearson"))
    
  }
  
  corScores <- c(corScores, cor(effect1, effect2, "everything", "pearson"))
  # calculate and reeturn pearson correlation
  return(corScores)
  
}

corDrug <- function (drug1, drug2) {
  
  corScores <- c()
  
  allDrugs <- getAllDrugs()
  
  dose <- seq (0, 1, 0.01)
  
  # dfm <- drugFeatureMatrix
  dfm <- getDrugFeatureMatrix()
  
  coefDRC1 <- dfm[match(drug1, allDrugs),]
  coefDRC2 <- dfm[match(drug2, allDrugs),]
  
  effect1 <- c()
  effect2 <- c()
  
  # for all cell lines
  for (i in seq(0,length(coefDRC1)/3 - 1)) {
    
    # interpolate data points from the curves
    newEffects1 <- fHill(dose, coefDRC1[i*3 + 1], coefDRC1[i*3 + 2], coefDRC1[i*3 + 3])
    newEffects2 <- fHill(dose, coefDRC2[i*3 + 1], coefDRC2[i*3 + 2], coefDRC2[i*3 + 3])
    effect1 <- c(effect1, newEffects1)
    effect2 <- c(effect2, newEffects2)
    
    corScores <- c(corScores, cor(newEffects1, newEffects2, "everything", "pearson"))
    
  }
  
  corScores <- c(corScores, cor(effect1, effect2, "everything", "pearson"))
  # calculate and return pearson correlation
  return(corScores)
  
}

# == COMPUTE SIMILARITY MATRICES ========================================================
# For drug: using mono therapy data and phy/chem data
# For cell line: using mono therapy data only (for now)

#Euclidean sum
norm_vec <- function(x){
  print(x)
  return sqrt(sum(x^2, na.rm=TRUE))
}
# call several sub functions to compute distance score between drug pairs for each cell line.  
computeDrugDistanceMatrix <- function(weights) {
  
  drugPropertiesWeight <- 0.2
  
  drugMonoWeight<- 1-drugPropertiesWeight
  
  druglist <-  read.csv(DRUG_PROPERTIES, head=TRUE,stringsAsFactors=FALSE)[1]
  drugMatrix <-data.frame(matrix(NA, nrow = nrow(druglist), ncol = nrow(druglist)))
  colnames(drugMatrix) <- druglist["ChallengeName"][,1]
  rownames(drugMatrix) <- druglist["ChallengeName"][,1]
  
  celllist <-  read.csv(CELL_PROPERTIES, head=TRUE,stringsAsFactors=FALSE)[1]
  cellMatrix <-data.frame(matrix(NA, nrow = nrow(celllist), ncol = nrow(celllist)))
  colnames(cellMatrix) <- celllist["Sanger.Name"][,1]
  rownames(cellMatrix) <- celllist["Sanger.Name"][,1]
  
  # Get the distance based off of physical / drug properties
    for(d1 in druglist[,1]){
      drugMatrix[d1,d1]<-1
      for(d2 in druglist[,1]){
        val <- norm_vec(distDrug(getDrugProperties(d1),getDrugProperties(d2)))
        drugMatrix[d1,d2] <- val
        drugMatrix[d2,d1] <- val
      }
    }
    
    # Normalize Drugmatrix 0 to 1.
    drugMatrix <- t(apply(drugMatrix, 1, function(x)(x-min(x))/(max(x)-min(x))))
    
    # Change distance to similarity.
    drugMatrix <- 1-drugMatrix
  
  print("Processing Drugs")
  for(d1 in druglist[,1]){
    drugMatrix[d1,d1]<-1
    for(d2 in druglist[,1]){
      print(c(d1, d2))
      drugMatrix[d1,d2]<-mean(corDrug(d1,d2), na.rm=TRUE)
      drugMatrix[d2,d1]<-drugMatrix[d1,d2]
    }
  }
  print("Processing Cells")
  for(c1 in cellist){
    cellMatrix[c1,c1]<-1
    for(c2 in cellist){
      print(c(c1, c2))
      cellMatrix[c1,c2]<- mean(corCellLine(c1,c2), na.rm=TRUE)
      cellMatrix[c2,c1]<-cellMatrix[c1,c2]
    }
  }
  
  write.table(drugMatrix, file = "DSM.csv")
  write.table(cellMatrix, file = "CSM.csv")
  # Manhattan distance idea: distance is the sum of distances of their corresponding components 
  # as if traveling on a grid (i.e. simply sum up the distances in each dimension)
  
  # Euclidean 
  # Straight line distance from datapoint
  # square root of sum of squares
  
  
}

# CELL LINE DEPENDENT FEATURES
# 2 main method of comparing correlation:
#    1. Looking at one cell line at a time 
#        i. curve fit mono therapy drug response to drug A and drug B
#        ii. interpolate n data points from each curve
#        iii. plot kill at concentrations (one drug on an axis) and compute correlation
#    2. Looking at all the cell lines together
#        Ideally: we want one value representing the effectiveness of a drug on a cell line
#                 i.e. want to combine IC50, Einf and H, but how to combine them in a meaningful way?
#        For now:
#        for drug_pair in all_drug_pairs[]:
#           for cell_line in all_cell_lines[]:
#               for property in (IC50, Einf, H):
#                   plot (observed value of property of drug_pair[0] on cell_line, observed value of property of drug_pair[1] on cell_line)


corDrugGivenCellLine <- function (drug1, drug2, cellLine) {
  
  if(drug1==drug2){
    return(1)
  }
  
  print (c(drug1,drug2,cellLine))
  
  dose <- seq (0, 1, 0.01)
  
  # stores into a data frame, mono therapy effects given a specific drug and cell line
  dfMonoEffect1 <- getMonoEffect(drug1, cellLine)
  dfMonoEffect2 <- getMonoEffect(drug2, cellLine)
  
  if(is.null(dfMonoEffect1[8:13]) || is.null(dfMonoEffect2[8:13])){
    return(0)
  }
  
  # get drug response curves 
  drc1 <- list()
  drc2 <- list()
  # normalize the concentrations
  drc1$conc <- c(0, 0.01, 0.03, 0.1, 0.3, 1)
  drc2$conc <- c(0, 0.01, 0.03, 0.1, 0.3, 1)
  # take the average of the effects 
  print(c(dfMonoEffect1[8:13], dfMonoEffect2[8:13]))
  
  drc1$dat <- colMeans(dfMonoEffect1[8:13])
  drc2$dat <- colMeans(dfMonoEffect2[8:13])
  
  # curve fit 
  coefDRC1 <- nlsDRC(drc1)
  coefDRC2 <- nlsDRC(drc2)
  
  # interpolate data points from the curves
  effect1 <- fHill(dose, coefDRC1[1], coefDRC1[2], coefDRC1[3])
  effect2 <- fHill(dose, coefDRC2[1], coefDRC2[2], coefDRC2[3])
  
  # calculate and return pearson correlation
  return(cor(effect1, effect2, "everything", "pearson"))
  
}

# compute and return correlation of two cell lines given a drug
corCellLineGivenDrug <- function (cellLine1, cellLine2, drug) {
  
  if(cellLine1==cellLine2){
    return(1)
  }
  
  dose <- seq (0, 1, 0.01)
  
  # stores into a data frame, mono therapy effects given a specific drug and cell line
  dfMonoEffect1 <- getMonoEffect(drug, cellLine1)
  dfMonoEffect2 <- getMonoEffect(drug, cellLine2)
  
  # get drug response curves 
  drc1 <- list()
  drc2 <- list()
  # normalize the concentrations
  drc1$conc <- c(0, 0.01, 0.03, 0.1, 0.3, 1)
  drc2$conc <- c(0, 0.01, 0.03, 0.1, 0.3, 1)
  # take the average of the effects 
  drc1$dat <- colMeans(dfMonoEffect1[8:13])
  drc2$dat <- colMeans(dfMonoEffect2[8:13])
  
  # curve fit 
  coefDRC1 <- nlsDRC(drc1)
  coefDRC2 <- nlsDRC(drc2)
  
  # interpolate data points from the curves
  effect1 <- fHill(dose, coefDRC1[1], coefDRC1[2], coefDRC1[3])
  effect2 <- fHill(dose, coefDRC2[1], coefDRC2[2], coefDRC2[3])
  
  # calculate and return pearson correlation
  return(cor(effect1, effect2, "everything", "pearson"))
  
}


#drugFeatureMatrix <- getDrugFeatureMatrix()
#cellFeatureMatrix <- getCellFeatureMatrix()

# Questions:
# 1. Observed drug effect as function of drug concentration (Chapter 2.2). When a = conc that gives Einf.... the IC50/a *or^H should be zero.... but how?!
# 2. Discrepancies between mono therapy observed effects.... which to use? Rule out outliners and take the mean?
