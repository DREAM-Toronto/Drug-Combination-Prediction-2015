PATHWAYS_FILE <- "c6.all.v5.0.symbols.gmt"
DRUGS_FILE <- "../Challenge Data/Drug Synergy Data/Drug_info_release.csv"
GEX_FILE <- "../Challenge Data/Sanger Molecular Data/gex.csv"

NCELLS <- 83 #number of cell lines in GEX_FILE (ncols - 1)

HIGH_GEX = 9   ###how to determine threshold values for high/low gene expression for cells?
LOW_GEX = 0.2

if (! require(magrittr, quietly=TRUE)) {
  # for plotting dose-response surfaces
  install.packages("magrittr",
                   repos="http://cran.us.r-project.org")
  library(magrittr)
}

# ================================================================================
readGEX <- function(gex_file = GEX_FILE) {
  gex <- read.csv(gex_file, header = TRUE,
                  colClasses = c("character", rep("numeric", NCELLS)),
                  check.names = FALSE, stringsAsFactors = FALSE)
  return(gex)
}

# ================================================================================
# get list of genes targeted by each drug
readDrugs <- function(drugs_file = DRUGS_FILE) {
  raw <- read.csv(drugs_file, header = TRUE, stringsAsFactors = FALSE)
  drugs <- raw[, 1]
  targets <- raw[, 2]
  result <- list()
  
  for(i in 1:nrow(raw)) {
    tmp <- strsplit(targets[i], ",")
    result <- c(result, tmp)
  }
  
  names(result) <- drugs
  return (result)
}

# ================================================================================
readPathways <- function(pathways_file = PATHWAYS_FILE) {
  raw <- read.csv("c6.all.v5.0.symbols.csv",
                  head=FALSE,
                  stringsAsFactors=FALSE)
}

# ================================================================================
# returns all genes in the same pathway as the drug target genes
# as a character vector with no duplicates
getDrugRelatedGenes <- function(drug_name, drugs, pathways) {
  targets <- drugs[[drug_name]]
  result <- c()
  for(i in 1:length(targets)) {
    result <- c(result, getRelatedGenes(targets[i], pathways))
  }
  return (unique(result))
}

# ================================================================================
# returns all genes in the same pathway as the highest expressed genes in the cell,
# another vector for those of the lowest expressed genes in the cell
getCellRelatedGenes <- function(cell_name, gex, pathways) {
  cell_gex <- gex[, cell_name]
  names(cell_gex) <- gex[, 1]
  max <- max(cell_gex)
  min <- min(cell_gex)
  
  high_genes <- c()
  low_genes <- c()
  
  for(i in 1:length(cell_gex)) {
    cell_gex[i] <- (cell_gex[i] - min) / ((max - min) / 10) #normalize values to 0-10
    if(cell_gex[i] >= HIGH_GEX) {
      high_genes <- c(high_genes, names(cell_gex)[i])
    } else if(cell_gex[i] <= LOW_GEX) {
      low_genes <- c(low_genes, names(cell_gex)[i])
    }
  }
  
  high_related <- c()
  low_related <- c()
  
  for(i in 1:length(high_genes)) {
    high_related <- c(high_related, getRelatedGenes(high_genes[i], pathways))
  }
  for(i in 1:length(low_genes)) {
    low_related <- c(low_related, getRelatedGenes(low_genes[i], pathways))
  }
  
  return(list(high = unique(high_related), low = unique(low_related)))
}

# ================================================================================
# returns all genes in the same pathways as the given gene
getRelatedGenes <- function(gene, pathways) {
  # TODO gene = is read from drug_info_release.csv = the drug target but can be "DNA", "proteosome", etc.??
  
  # genes can be listed with '*' denoting any character if there is uncertainty as to the specific gene i.e. AKT* can be AKT1, AKT2...
  if(length(grep("\\*", gene) > 0)) {
    gene <- gsub("\\*", "(.)*", gene)
  }
  
  containedGenes <- c()
  nRow <- nrow(pathways)
  nCol <- ncol(pathways)
  for(i in 1:nRow){
    currRow <- c(pathways[i,3:nCol], strsplit(pathways[i, 1], '_')[[1]][1]) #adds gene in col1
    # TODO sometimes in pathways file, ex: gene entered from col1 is AKT instead of AKT1 or AKT*, and is not matched if gene = AKT1
    if(length(grep(gene, currRow) > 0)){
      containedGenes <- union(containedGenes, currRow)
    }
  }
  return (containedGenes)
}

# ================================================================================
# calculate intersect / union of genes1 and genes2
calculate_overlap <- function(genes1, genes2) {
  return (length(intersect(genes1, genes2)) / length(union(genes1, genes2)))
}

# ================================================================================

gex <- readGEX()
drugs <- readDrugs()
pathways <- readPathways()
therapy <- read.csv("../Challenge Data/Drug Synergy Data/ch2_leaderBoard_monoTherapy.csv", header = TRUE, stringsAsFactors = FALSE)
pathway_features <- matrix(nrow = nrow(therapy), ncol = 8,
                           dimnames = list(c(), c("CELL_LINE", "COMPOUND_A", "COMPOUND_B",
                                             "HIGH, COMPOUND_A", "LOW, COMPOUND_A",
                                             "HIGH, COMPOUND_B", "LOW, COMPOUND_B",
                                             "COMPOUND_A, COMPOUND_B")))

for(i in 1:nrow(therapy)) {
  cell <- therapy[i, 1]
  drug1 <- therapy[i, 2]
  drug2 <- therapy[i, 3]
  
  cell_changed <- TRUE; drug1_changed <- TRUE; drug2_changed <- TRUE;
  if(i > 1) {
    if(cell == therapy[i - 1, 1]) cell_changed = FALSE
    if(drug1 == therapy[i - 1, 2]) drug1_changed = FALSE
    if(drug2 == therapy[i - 1, 3]) drug2_changed = FALSE
  }
  
  print(sprintf("Making features for %d: %s x %s x %s", i, cell, drug1, drug2))
  
  if(cell_changed) cell_pathway_genes <- getCellRelatedGenes(cell, gex, pathways)
  if(drug1_changed) drug1_pathway_genes <- getDrugRelatedGenes(drug1, drugs, pathways)
  if(drug2_changed) drug2_pathway_genes <- getDrugRelatedGenes(drug2, drugs, pathways)
  
  pathway_features[i, 1] = cell
  pathway_features[i, 2] = drug1
  pathway_features[i, 3] = drug2
  pathway_features[i, 4] = calculate_overlap(cell_pathway_genes$high, drug1_pathway_genes)
  pathway_features[i, 5] = calculate_overlap(cell_pathway_genes$low , drug1_pathway_genes)
  pathway_features[i, 6] = calculate_overlap(cell_pathway_genes$high, drug2_pathway_genes)
  pathway_features[i, 7] = calculate_overlap(cell_pathway_genes$low , drug2_pathway_genes)
  pathway_features[i, 8] = calculate_overlap(drug1_pathway_genes    , drug2_pathway_genes)
}

#END
