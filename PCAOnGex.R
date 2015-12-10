setwd("~/Dream_Challenge/Drug-Combination-Prediction-2015")
gex_file <- "../Challenge\ Data/Sanger\ Molecular\ Data/gex.csv"
raw_cellNames <- read.csv(gex_file, header = FALSE, stringsAsFactors = FALSE,nrows = 1)
raw_cellNames <- as.character(raw_cellNames)
raw_cellNames <- raw_cellNames[2:length(raw_cellNames)]
raw <- read.csv(gex_file, header = TRUE, stringsAsFactors = FALSE)
rownames(raw) <- raw[,1]
raw[,1] <- NULL
colnames(raw) <- raw_cellNames

raw_pca <- prcomp(raw,center = TRUE,scale. = TRUE)
write.csv(raw_pca$rotation,"gex_pca.csv")