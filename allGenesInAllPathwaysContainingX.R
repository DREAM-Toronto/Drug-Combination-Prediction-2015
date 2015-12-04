
allGenesInAllPathwaysContainingX = function(gene){
  raw <- read.csv("c6.all.v5.0.symbols.csv",
                  head=FALSE,
                  stringsAsFactors=FALSE)
  containedGenes <- c()
  nRow <- nrow(raw)
  nCol <- ncol(raw)
  for(i in 1:nRow){
    currRow <- raw[i,3:nCol]
    if(gene %in% currRow == TRUE){
      containedGenes <- union(containedGenes, raw[i,3:nCol])
    }
  }
  return (containedGenes)
}


