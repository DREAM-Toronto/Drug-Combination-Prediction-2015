# TO-DO's:
# - Fill in missing information in Drug_info_release.csv
# - Implement cell line response comparison (AUC instead of IC50,H,Einf?)
# - Impelement a main function: computeDrugSimilarityMatrix (see below)
# - Try different weights for different features (more weights on cell line dependent features?)
# - Allow choice of different similarity metrics for cell line dependent features

MONO_FILE <- "../Challenge Data/Drug Synergy Data/Drug_info_release.csv"
DRUG_PROPERTIES <- "../Challenge Data/Drug Synergy Data/Drug_info_release.csv"
CELL_KILL_DIR <- "../Challenge Data/Drug Synergy Data/Raw Data/Raw_Data_csv/ch1_ch2_monoTherapy/"
CELL_PROPERTIES <- "../Challenge Data/Sanger Molecular Data/cell_info.csv"

getMonoEffect <- function(drug,cellLine) {
  CELL_KILL_DIR <- "../Challenge Data/Drug Synergy Data/Raw Data/Raw_Data_csv/ch1_ch2_monoTherapy/"
  
  regEx <- paste(".*(",drug,"[.]).*(",cellLine,")",sep="")
  filenameList <- list.files(CELL_KILL_DIR,pattern = regEx)
  
  allMonoEffects <- c()
  
  for(fn in filenameList){
    raw <- read.csv(paste(CELL_KILL_DIR, fn, sep=""),
                    head=FALSE,
                    stringsAsFactors=FALSE)
    
    filename <- fn
    conc <- "null"
    effectiveness <- "null"
    monoEffect <- list()
    
    if(raw[9,2]==drug){
      # Agent1
      conc <- raw[2:7,1]
      effectiveness <- raw[2:7,2]
    }else{
      # Agent2
      conc <- raw[1,2:7]
      effectiveness <- raw[2,2:7]
    }
    monoEffect <- list()
    monoEffect$FN <- filename
    monoEffect$Conc <- conc
    monoEffect$Eff <- effectiveness
    
    allMonoEffects <- c(allMonoEffects,monoEffect)
  }
  
  return(allMonoEffects)
}


# == PARSE DRUG INFO CSV ========================================================

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

distIC50 <- function(IC501, IC502) {
  abs(IC501 - IC502)  
}

distEinf <- function(Einf1, Einf2) {
  abs(Einf1 - Einf2)  
}

distH <- function(H1, H2) {
  abs(H1 - H2)  
}

# compute the distance between two drugs
distDrug <- function(drug1, drug2) {

# cell line dependent features 
# will separate cell line features from the other features
#   for (i in 1:85) {
#     cellLineResponse <- list(
#       cellLineResponse, list(
#       distIC50(drug1[(i-1)*3 + 1], drug2[(i-1)*3 + 1]),
#       distEinf(drug1[(i-1)*3 + 2], drug2[(i-1)*3 + 2]),
#       distH(drug1[(i-1)*3 + 3], drug2[(i-1)*3 + 3]))
#       )    
#   }
  # These distances need a weight...
  distance <- c(
    distTarget(drug1[[1]],drug2[[1]]),
    distHBA(drug1[[2]],drug2[[2]]),
    distCLogP(drug1[[3]],drug2[[3]]),
    distHBD(drug1[[4]],drug2[[4]]),
    distMW(drug1[[5]],drug2[[5]])
#    cellLineResponse
  )
  return(distance)
}

#Euclidean sum
norm_vec <- function(x) sqrt(sum(x^2, na.rm=TRUE))

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
#   for(d1 in druglist[,1]){
#     drugMatrix[d1,d1]<-1
#     for(d2 in druglist[,1]){
#       val <- norm_vec(distDrug(getDrugProperties(d1),getDrugProperties(d2)))
#       drugMatrix[d1,d2] <- val
#       drugMatrix[d2,d1] <- val
#     }
#   }
#   
#   # Normalize Drugmatrix 0 to 1.
#   drugMatrix <- t(apply(drugMatrix, 1, function(x)(x-min(x))/(max(x)-min(x))))
#   
#   # Change distance to similarity.
#   drugMatrix <- 1-drugMatrix
  
  for(d1 in druglist[,1]){
    drugMatrix[d1,d1]<-1
    for(d2 in druglist[,1]){
      values <- c()
      for(cl in celllist[,1]){
        values <- c(values,corDrugGivenCellLine(d1,d2,cl))
      }
      drugMatrix[d1,d2]<-drugPropertiesWeight*drugMatrix[d1,d2] +  drugMonoWeight*mean(values)
      drugMatrix[d2,d1]<-drugMatrix[d1,d2]
    }
  }
  
  for(c1 in cellist){
    cellMatrix[c1,c1]<-1
    for(c2 in cellist){
      values <- c()
      for(d in druglist){
        values <- c(values)
      }
      cellMatrix[c1,c2]<-mean(values)
      cellMatrix[c2,c1]<-mean(values)
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

# Questions:
# 1. Observed drug effect as function of drug concentration (Chapter 2.2). When a = conc that gives Einf.... the IC50/a *or^H should be zero.... but how?!
# 2. Discrepancies between mono therapy observed effects.... which to use? Rule out outliners and take the mean?