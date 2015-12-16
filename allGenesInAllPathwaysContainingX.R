
allGenesInAllPathwaysContainingX = function(gene){
  raw <- read.csv("c6.all.v5.0.symbols.csv",
                  head = TRUE,
                  stringsAsFactors = FALSE)
  containedGenes <- c()
  nRow <- nrow(raw)
  nCol <- ncol(raw)
  for(i in 1:nRow){
    currRow <- unlist(raw[i,3:nCol]) 
    currRow <- currRow[currRow!=""] # ignore empty strings
    if(gene %in% currRow == TRUE){
      containedGenes <- union(containedGenes, currRow)
    }
  }
  return (containedGenes)
}


