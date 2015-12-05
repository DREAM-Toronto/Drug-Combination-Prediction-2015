setwd(DREAMDIR)

PATHWAYS_FILE <- "../Challenge Data/External/c6.all.v5.0.symbols.gmt"
DRUGS_FILE <- "../Challenge Data/Drug Synergy Data/Drug_info_release.csv"
GEX_FILE <- "../Challenge Data/Sanger Molecular Data/gex.csv"

HIGH_GEX = 9
LOW_GEX = 0.2

if (! require(magrittr, quietly=TRUE)) {
  # for plotting dose-response surfaces
  install.packages("magrittr",
                   repos="http://cran.us.r-project.org")
  library(magrittr)
}

# ================================================================================
readGEX <- function(gex_file = GEX_FILE) {
  raw <- read.csv(gex_file, header = FALSE, stringsAsFactors = FALSE)
  cell_lines <- raw[1,2:84]
  genes <- raw[2:17420,1] 
  
  tmp <- raw[2:17420, 2:84]%>%
    unlist %>%
    as.numeric
  gex <- matrix(tmp, nrow = 17419, ncol = 83)
  dimnames(gex) <- list(genes, cell_lines)
  return(gex)
}

# ================================================================================
readDrugs <- function(drugs_file = DRUGS_FILE) {
  
  
}

# ================================================================================
readPathways <- function(pathways_file = PATHWAYS_FILE) {
  
   
}

# ================================================================================
# returns all genes in the same pathway as the drug target genes
getDrugRelatedGenes <- function(drug_name, drugs, pathways) {
  
  
}

# ================================================================================
# returns all genes in the same pathway as the highest expressed genes in the cell,
# another vector for those of the lowest expressed genes in the cell
getCellRelatedGenes <- function(cell_name, gex, pathways) {
  cell_gex <- gex[, cell_name]
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
    high_related <- c(high_related, getRelatedGenes(high_genes[i]), pathways)
  }
  for(i in 1:length(low_genes)) {
    low_related <- c(low_related, getRelatedGenes(low_genes[i]), pathways)
  }
  
  return(list(high = unique(high_related), low = unique(low_related))) ######get rid of duplicates??
}

# ================================================================================
# returns all genes in the same pathways as the given gene
getRelatedGenes <- function(gene, pathways) {
  #implement...
  
  
  return(gene)
  
}

# ================================================================================
# calculate intersect / union of genes1 and genes2
calculate_overlap <- function(genes1, genes2) {
  
}

# ================================================================================
# returns a matrix of #drug combos x 3
# columns are % overlap of related genes of
# 1. drug1 targets  & drug2 targets
# 2. highest and 3. lowest expressed genes in cell & drug targets
makePathwayFeatures <- function(combo_therapy_file, gex_file = GEX_FILE, drugs_file = DRUGS_FILE, pathways_file = PATHWAYS_FILE) {
  #implement...
  
  gex <- readGEX(gex_file)
  drugs <- readDrugs(drugs_file)
  pathways <- readPathways(pathways_file)
  
  
  
}

# END
