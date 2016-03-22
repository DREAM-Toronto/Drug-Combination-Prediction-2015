# makeMSig_GEXFeatures.R
#
# Purpose: Compile features for drug x drug and drug x drug x cell
#          combinations by using MSigDB gege sets, Sanger
#          expression data and a STRING functional relationship
#          network.
#
# Version: 1.0
#
# Date:    March 21 2016
# Author:  Boris and DREAM team UofT
#
# V 1.0    Implemented Expression weigthing - but not in time for the
#          submission.
# V 0.2    Implemented Fuzzy Jaccard Distance
# V 0.1    Jaccard Distance only
#
# TODO:
#          Compute drug x drug x cell scores and save the hash
#          Explore correlation with synergy scores
# ====================================================================


setwd("/Users/steipe/Documents/00.3.REFERENCE_AND_SUPPORT/DREAM 2015/work")

library(igraph)
library(hash)
if (!require(biomaRt)) {
  source("http://bioconductor.org/biocLite.R")
  biocLite("biomaRt")
  library("biomaRt")
}

# === FILE NAMES AND LOCATIONS =======================================

# Gene expression data provided by Challenge authors
GEXFILE <- "../Challenge Data/Sanger Molecular Data/gex.csv"

# Drug information data provided by Challenge authors, modified by
# us to augment the list of drug targets
DRUGFILE <- "../Challenge Data/Drug Synergy Data/Drug_info_release.mod.csv"

# "Pathways" downloaded from MSigDB. C6 is "cancer pathways"
MSIGFILE <- "~/Documents/00.3.REFERENCE_AND_SUPPORT/DREAM 2015/work/c6.all.v5.0.symbols.csv"

# Drug dose/response/synergy data provided by Challenge authors
DRSFILE <- "../Challenge Data/Drug Synergy Data/ch1_train_combination_and_monoTherapy.csv"

# If you don't have the STRING file for humans (9606), download a version
# from http://string.embl.de/newstring_cgi/show_download_page.pl
# ... it's about 500MB
STRINGFILE <- "/Users/steipe/Documents/07.TEACHING/50.8-BCB420-JTB2020 2016/BCB420/data/STRING/9606.protein.links.detailed.v10.txt"

# === READ AND PROCESS FILES =========================================

# --------  drugs
tmp <- read.csv(DRUGFILE,
                header = TRUE,
                stringsAsFactors = FALSE)
drugs <- list()

for(i in 1:nrow(tmp)) {
    s <- unlist(strsplit(tmp$Target.Official.Symbol[i], "\\s*,\\s*"))
    drugs[[tmp$ChallengeName[i]]] <- s
}
# length(drugs)
# drugs is a list of length 119. Each list element is named as a drug
# and has the value of a character vector containing one or
# more target gene symbols.

# --------  geneSets
tmp <- read.csv(MSIGFILE,
                head=FALSE,
                stringsAsFactors=FALSE)
C6set <- list()
for(i in 1:nrow(tmp)) {
  s <- unname(unlist(tmp[i, tmp[i,] != ""]))
  C6set[[i]] <- s
}
# length(C6set)
# C6set is a list of length 189. Each list element
# has the value of a character vector containing the
# target gene symbols for a "cancer gene pathway" from
# the MSigDB C6 collection.


# --------  dose/response/synergy data
drugEffectData <- read.csv(DRSFILE,
                           stringsAsFactors=FALSE)
drugEffectData <- drugEffectData[drugEffectData$QA == 1, ] # Exclude all combinations with poor QA scores



# ====================================================================
# PART 1:  JACCARD DISTANCE
# ====================================================================


# ======  HELPER FUNCTIONS AND DEFINITIONS ===========================

makeGrepPattern <- function(s) {
  # Given a symbol with possible wildcard expansion
  # make a regex pattern that will identify it in a symbol
  #
  if (length(grep('\\*', s) > 0)) {
    return(sprintf("^%s$", gsub('\\*', '.*', s))) # eg. "^AKT.*"
  } else {
    return(sprintf("^%s$", s)) # eg. "^ADAM17$"
  }
}


# Calculating set overlaps is slow, but repetitive. I have achieved
# dramatic speedup from hashing all results.

PIhash <- hash()
getPathwayIndices <- function(targets, set, reset=FALSE) {
  # return the indices of all elements of "set" that
  # contain an element of "targets" or its expansion
  if (reset) {
    PIhash <<- hash()
    return()
  }
  key <- paste(targets, collapse = ".")
  paths <- PIhash[[key]]
  if ( is.null(paths)) {
    paths = numeric()
    for (i in 1:length(targets)) {
      pattern <- makeGrepPattern(targets[i])
      x <- lapply(set, function(x) as.logical(length(grep(pattern, x) > 0))  )
      x <- which(unlist(x))

      if (length(x) > 0) {
        paths <- c(paths, x)
      }
    }
    paths <- unique(paths)
    PIhash[[key]] <<- paths
  }
  return(paths)
}


PGhash <- hash()
getPathwayGenes <- function(idx, set, reset=FALSE) {
  # return the genes of all elements of "set" that
  # are referenced in "idx"
  if (reset) {
    PGhash <<- hash()
    return()
  }
  key <- sprintf("X%s",paste(idx, collapse = "."))
  genes <- PGhash[[key]]
  if ( is.null(genes)) {
    genes = character()
    if (length(idx) > 0) {
      for (i in 1:length(idx)) {
        genes <- c(genes, set[[idx[i]]])
      }
      genes <- unique(genes)
    }
    PGhash[[key]] <<- genes
  }
  return(genes)
}


dJaccard <- function(A, B) {
  # Jaccard Distance is 1-Jaccard Index
  return(1 - length(intersect(A, B)) / length(union(A, B)))
}


JDhash <- hash()
JD <- function(A, B, drugs, pathways, reset=FALSE) {
  # This calculates the Jaccard Distance for drug targets pathways.
  # A: string: compound A
  # B: string: compound B
  # drugs: a list of drug targets
  # pathways: a list of pathways that contain drug targets.
  #
  # This function maintains its results in a hash.
  # Call the function with JD(reset=TRUE) to reset the hash.
  #
  # The Jaccard Distance is 1 - intersect(X, Y) / union(X, Y) where
  # X and Y are the unions of pathway target genes for compound A
  # and B respectively.
  #
  # If a target is given as XYZ* it is matched with all symbols that
  # begin with XYZ.
  #
  if (reset) {
    JDhash <<- hash()
    return()
  }
  key <- sprintf("%s.%s", A, B)
  dJ <- JDhash[[key]]
  if (is.null(dJ)) {  # not present in hash ...
    aIdx <- getPathwayIndices(drugs[[A]], pathways)
    bIdx <- getPathwayIndices(drugs[[B]], pathways)
    if (length(aIdx) == 0 | length(bIdx) == 0 ) {
      dJ <- 1
    } else {
      aSet <- getPathwayGenes(aIdx,  pathways)
      bSet <- getPathwayGenes(bIdx,  pathways)
      dJ <- dJaccard(aSet, bSet)
    }
    JDhash[[key]] <<- dJ      # store in environment
  }
  return(dJ)
}


# Plot the Jaccard Distance for all drugs
# (Analytics only)
# N <- length(drugs)
# cat("\n")
# dMat <- matrix(numeric(N * N), nrow=N)
# for (i in 1:N) {
#   cat("=")
#   for (j in i:N) {
#     dMat[i, j] <- dMat[j, i] <- JD(names(drugs)[i], names(drugs)[j], drugs, C6set)
#   }
#   if (!(i %% 20)) {cat("\n")}
# }
# cat("\n")
# image(dMat, col = colorRampPalette(c("#000000", "#5588FF", "#FFFFEE"))(12))


# ======  CALCULATING FEATURES =======================================

# The Jaccard distance can be used to calculate drug x drug similarity.
# use as:

# value <- JD(A, B, drugs, C6set)

# where:
# A: is COMPOUND_A
# B: is COMPOUND_B
# and drugs and C6set have been defined above

# Example:
JD(drugEffectData$COMPOUND_A[123], drugEffectData$COMPOUND_B[123], drugs, C6set) # 0.9906542


# ====================================================================
# PART 2:  FUZZY JACCARD DISTANCE
# ====================================================================

# Fuzzy Jaccard distances are calculated by adding scaled scores for
# overlapping neighborhoods. Neighborhoods are taken from the STRING
# graoh of functional relationships.

# === READ AND PROCESS FILES =========================================

# --------  STRING graph
#

tmp <- read.delim(STRINGFILE, header=TRUE, sep=" ", stringsAsFactors=FALSE)
# nrow(tmp)  # 8,548,002
tmp$protein1 <- substr(tmp$protein1, 6, 20)  # drop the "9606." prefix
tmp$protein2 <- substr(tmp$protein2, 6, 20)

ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

allGenes <- unique(tmp$protein1)
esMap <- getBM(filters = "ensembl_peptide_id",    # about 20 sec.
               attributes = c("ensembl_peptide_id",
                              "hgnc_symbol"),
               values = allGenes,
               mart = ensembl)

# length(allGenes)  # 19,247
# nrow(esMap)       # 18,168
colnames(esMap) <- c("ens", "sym")
rownames(esMap) <- esMap$ens
# head(esMap)

# substitute gene symbols for ENS IDs
tmp$protein1 <- esMap[tmp$protein1, "sym"]   # about 1 min. each
tmp$protein2 <- esMap[tmp$protein2, "sym"]   #

# drop all rows in which either protein is mapped to "" or NA
#
tmp <- tmp[(tmp$protein1 != "" &
            tmp$protein2 != "" &
            !(is.na(tmp$protein1)) &
            !(is.na(tmp$protein2))), ]
# nrow(tmp)  # 7,680,670

thrsh <- 900  # Empiricaly determined by considering overlap fraction
              # of pathway genes and STRING graph, confidence, and
              # biological parameters. At 900, we have an average node
              # degree of ~ 10 - this rises to 20 for e.g. 800 and
              # that probably implies that we have 50% false positive
              # edges...
              #
              # Perhaps an improvement could be obtained by further
              # connecting all remaining pathway members along their
              # highest-confidence path? To test on a rainy day ...

STRINGedges <- tmp[tmp$combined_score > thrsh, ]
# nrow(STRINGedges)  # 221,880
# nodes <- unique(c(STRINGedges$protein1, STRINGedges$protein2))
# length(nodes)  # 10,440

# Do these contain all of our drugTargets?
# MSigNodes <- unique(unlist(C6set, use.names=FALSE))
# length(MSigNodes)  # 11,250
# length(intersect(MSigNodes, nodes))  # 6,641
# length(intersect(MSigNodes, nodes)) / length(nodes)  # 63.6 % of nodes in the graph
# nrow(STRINGedges) / length(nodes)  # Average degree is 0.5 of 21.3 ... that's about right
rm(tmp)

gSTR <- graph_from_data_frame(STRINGedges[ , c(1, 2, 10)],
                              directed = FALSE) # create the igraph object

# Here we keep only the combined_score as edge attribute, we don't
# actually use the score however, but treat gSTR as an unweighted,
#  undirected graph.

gSTRnodes <- unique(c(STRINGedges$protein1, STRINGedges$protein1))

# save(gSTR, gSTRnodes, file = "STRINGgraph.Rdata")
# load("STRINGgraph.Rdata")


# ======  HELPER FUNCTIONS AND DEFINITIONS ===========================


# scalF provides relative scales for fuzzy set overlap regions so
# that relative weights between regions have a constant factor and
# the sum over all scaling factors is 1. We use this to globally
# scale the Jaccard distances for fuzzy regions
scalF <- function(q) {
  # Scaling of overlap regions for fuzzy set overlap
  p  <- numeric()
  pp <- numeric()

  p[1] <- 1
  p[2] <- p[1] * q
  p[3] <- p[2] * q

  pp[1] <- p[1] * p[1]
  pp[2] <- p[1] * p[2]
  pp[3] <- p[2] * p[1]
  pp[4] <- p[1] * p[3]
  pp[5] <- p[3] * p[1]
  pp[6] <- p[2] * p[2]
  pp[7] <- p[2] * p[3]
  pp[8] <- p[3] * p[2]
  pp[9] <- p[3] * p[3]

  pp <- pp/sum(pp)
  p <- sqrt(pp[c(1, 6, 9)])

 # cat(sprintf(" a:%1.3f                                      b:%1.3f                    c:%1.3f \n",
 #              p[1], p[2], p[3]))
 # cat(sprintf("aa:%1.3f ab:%1.3f ba:%1.3f ac:%1.3f ca:%1.3f bb:%1.3f bc:%1.3f cb:%1.3f cc:%1.3f \n\n",
 #             pp[1], pp[2], pp[3], pp[4], pp[5], pp[6], pp[7], pp[8], pp[9]))
  return(pp)
}

FSCALE <- scalF(1/5) # 0.650 (0.130 0.130) (0.026 0.026 0.026) (0.005 0.005) 0.001




# --------  analytics: running times

# # Average time for a node distance calculation:
# N <- 1000
# d <- numeric(N)
# ptm <- proc.time() # Start the stopwatch...
# for (i in 1:N) {
#   d[i] <- distances(gSTR, sample(V(gSTR), 1), sample(V(gSTR), N))
# }
# proc.time() - ptm  # How long did we take?
#                    # user  system elapsed
#                    # 17.032   0.700  17.818
# hist(d)
#
# # Average time for a node neighborhod calculation:
# N <- 1000
# d <- numeric(N)
# ptm <- proc.time() # Start the stopwatch...
# for (i in 1:N) {
#   d[i] <- length(unlist(ego(gSTR, 2, sample(V(gSTR), 1))))
# }
# proc.time() - ptm  # How long did we take?
#                    # user  system elapsed
#                    # 5.334   0.254   5.586
# hist(d)
#
# # Both algorithms are quite fast, but the processing of neighborhoods
# # is understandably quite a bit faster than pairwaise distance
# # calculations.


# === EXPAND PATHWAYS WITH NEIGHBORHOODS =============================


# Expand sets with neighborhoods

# reset hashes
getPathwayGenes(reset = TRUE)
getPathwayIndices(reset = TRUE)

# Make a copy of pathways that contains only genes that are also in
# the STRING graph

C6_STRset <- C6set

for (i in 1:length(C6_STRset)) {
  C6_STRset[[i]] <- C6_STRset[[i]][C6_STRset[[i]] %in% gSTRnodes]
}
pathways <- C6_STRset


# Create the list drugN that contains the pathway genes for all
# drug targets, and their 1- and 2- adjacent neighbourhoods

drugN <- list()

cat("\n")

for (i in 1:length(drugs)){
  cat("=")
  thisName <- names(drugs)[i]

  drugN[[thisName]]$N0 <- getPathwayGenes(getPathwayIndices(drugs[[thisName]], pathways),  pathways)

  if (length(drugN[[thisName]]$N0) == 0) {
    # Not in pathway, try use the targets themselves  by
    # finding them in gSTR nodes
    genes = character()
    for (j in 1:length(drugs[[thisName]])) {
      pattern <- makeGrepPattern(drugs[[thisName]][j])
      genes <- c(genes, gSTRnodes[grep(pattern, gSTRnodes)])
    }
    drugN[[thisName]]$N0 <- unique(genes)
  }

  # N1 neighbourhood
  x <- adjacent_vertices(gSTR, drugN[[thisName]]$N0)
  y <- names(V(gSTR)[unique(unlist(x, use.names=FALSE))])
  drugN[[thisName]]$N1 <- setdiff(y, drugN[[thisName]]$N0)

  # N2 neighbourhood
  x <- adjacent_vertices(gSTR, drugN[[thisName]]$N1)
  y <- names(V(gSTR)[unique(unlist(x, use.names=FALSE))])
  drugN[[thisName]]$N2 <- setdiff(y, union(drugN[[thisName]]$N0, drugN[[thisName]]$N1))

  if(! i %% 20) { cat("\n") }
}

cat("\n")

# save(drugN, file = "drugN_C6.Rdata")
# load("drugN_C6.Rdata")

# check ...
# for (i in 1:length(drugN)) {
#   n <- names(drugs)[i]
#   cat(sprintf("%s\t N0:%d\t N1:%d\t N2:%d\n",
#               n,
#               length(drugN[[n]]$N0),
#               length(drugN[[n]]$N1),
#               length(drugN[[n]]$N2)  ))
# }

# Calculate distance

FJDhash <- hash()
FJD <- function(A, B, drugs, pathN, scal = FSCALE, reset=FALSE) {
  # This calculates a Fuzzy Jaccard Distance for drug target pathway
  # neighborhoods.
  # A: string: compound A
  # B: string: compound B
  # drugs: a list of drug targets
  # pathN: a list of pathways that contain drug targets and their
  #        neighborhoods.
  #
  # This function maintains its results in a hash.
  # Call the function with FJD(reset=TRUE) to reset the hash.
  #
  # The Jaccard Distance is 1 - intersect(X, Y) / union(X, Y) where
  # X and Y are the unions of pathway target genes for compound A
  # and B respectively.
  # The Fuzzy Jaccard Distance uses weights to weight overlaps in
  # the N1 and N2 neighborhoods.
  #
  #
  if (reset) {
    FJDhash <<- hash()
    return()
  }
  key <- sprintf("%s.%s", A, B)
  dFJ <- FJDhash[[key]]
  if (is.null(dFJ)) {  # not present in hash ...
    a <- pathN[[A]]
    b <- pathN[[B]]
    if (length(a$N0) == 0 | length(b$N0) == 0) {
      dFJ <- 1
    } else {
      # calculate dFJ
      dFJ <-       scal[1] * dJaccard(a$N0, b$N0)
      dFJ <- dFJ + scal[2] * dJaccard(a$N0, b$N1)
      dFJ <- dFJ + scal[3] * dJaccard(a$N1, b$N0)
      dFJ <- dFJ + scal[4] * dJaccard(a$N0, b$N2)
      dFJ <- dFJ + scal[5] * dJaccard(a$N2, b$N0)
      dFJ <- dFJ + scal[6] * dJaccard(a$N1, b$N1)
      dFJ <- dFJ + scal[7] * dJaccard(a$N1, b$N2)
      dFJ <- dFJ + scal[8] * dJaccard(a$N2, b$N1)
      dFJ <- dFJ + scal[9] * dJaccard(a$N2, b$N2)
    }
    FJDhash[[key]] <<- dFJ      # store in hash
  }
  return(dFJ)
}


# Plot the Fuzzy Jaccard Distance for all drugs
# (Analytics only)
# N <- length(drugs)
# dFMat <- matrix(numeric(N * N), nrow=N)
# cat("\n")
# for (i in 1:N) {
#   cat("=")
#   for (j in i:N) {
#     dFMat[i, j] <- dFMat[j, i] <- FJD(names(drugs)[i], names(drugs)[j], drugs, drugN)
#   }
#   if (! (i %% 20)) { cat("\n")}
# }
# cat("\n")
# image(dFMat, col = colorRampPalette(c("#000000", "#5588FF", "#FFFFEE"))(12))
# image(dMat - dFMat, col = colorRampPalette(c("#FF0000", "#000000", "#00FF00"))(12))
# plot(dMat, dFMat, xlim=c(0,1), ylim=c(0,1))
#
#
# save(FJDhash, file = "FJDhash.Rdata")
# load("FJDhash.Rdata")
# save(drugs, file = "drugs.Rdata")
# load("drugs.Rdata")


# ======  CALCULATING FEATURES =======================================

# The Fuzzy Jaccard distance can be used just like the Jaccard Distance
# to calculate drug x drug similarity. It is somehwat smoother than
# the former.
# Use it as:

# value <- FJD(A, B, drugs, drugN)

# where:
# A: is COMPOUND_A
# B: is COMPOUND_B
# drugs is the list of drug targets
# drugN is the list of core pathway genes and their N-1 and N-2 neighborhoods

# Example:
#
FJD(drugEffectData$COMPOUND_A[123], drugEffectData$COMPOUND_B[123], drugs, drugN) # 0.9730703


# ====================================================================
# PART 3:  ADDING EXPRESSION WEIGHTS
# ====================================================================

# === READ AND PROCESS FILES =========================================

# --------  expression data
NCELLS <- 83
tmp <- read.csv(GEXFILE,
                header = TRUE,
                colClasses = c("character", rep("numeric", NCELLS)),
                check.names = FALSE, stringsAsFactors = FALSE)

gex <- as.matrix(tmp[ , -1], byrow=TRUE)
rownames(gex) <- tmp[ ,1]
# gex is a numeric matrix of 17,419 rows and 83 columns. Each value
# is the expression value of the gene named in the row in the cell
# line named in the column. Data description is here:
# https://www.synapse.org/#!Synapse:syn4231880/wiki/235651
# Values are Robust Multi-array Average (RMA) normalised with the
# R-package 'affy': raw intensity values are background corrected,
# log2 transformed and then quantile normalized.

# boxplot(gex[,1:83]) # distributions across cell-lines are very well
                      # comparable.

# To weight gene- importance, we will transform the RMA values to
# log-ratios of the median expression level

# ptm <- proc.time() # Start the stopwatch...
 logXq <- t(apply(gex, 1, function(x) {x - median(x)}))
# proc.time() - ptm  # How long did we take?
#                    #   user  system elapsed
#                    # 30.739   0.073  30.808

# # Some analytics ...
# boxplot(logXq[,1:83])  # Nicely distributed around 0 across columns
#
# # Count number of genes that change by less than 2 fold.
# baseX <- numeric()
# for (i in 1: ncol(logXq)) {
#   baseX[i] <- sum(abs(logXq[,i]) < 1)
# }
#
# nrow(logXq)
# plot(sort(baseX))
# nrow(logXq)    #   17,419
# summary(baseX) #    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#                #   14730   15440   15580   15570   15730   16120
# # The number of differentialy expressed genes varies by about 10%
# # across cell lines.
# #

# Selecting:
# We only needs expression values for genes that are actually in
# our pathway neighborhoods. We flatten our drugN list, unique() it
# and remove all genes that are not in that list.

pathGenes <- unique(unlist(drugN, use.names=FALSE))
selPathInX <- pathGenes %in% rownames(logXq)
sum(!selPathInX)  # we don't have expression values for 821 pathway genes.

selXInPath <- rownames(logXq) %in% pathGenes
sum(!selXInPath)  # 7,886 expression values are not in our pathways.

# First we drop the expression values that are not in the pathway.
tmp <- logXq[selXInPath, ]
nrow(tmp) # 9,533

# Then we add rows of zeros for all pathway genes that we are missing
n <- sum(!selPathInX)
x <- matrix(numeric(n * ncol(tmp)), nrow=n, ncol=ncol(tmp))
rownames(x) <- pathGenes[!selPathInX]
# head(x)

# Finally we add the block of zeros to the matrix.
tmp <- rbind(tmp, x)

# Now the rows of tmp have a one-to-one correspondence to our
# pathway genes.
# sum(! pathGenes %in% rownames(tmp))  # 0
# sum(! rownames(tmp) %in% pathGenes)  # 0

# Reweighting:
# We would like to replace the length() of union and intersection by
# a measure that weights differentially expressed genes more strongly.
# We scale each column: take the absolutes of the log(ratios) and
# rescale them so the column sums is the same as the number of rows.
# This means on average, each gene still contributes 1 to the
# cardinality of a set, but the actual values differ.
#

weightX <- apply(tmp, 2, function(x) {(abs(x) / sum(abs(x))) * length(x) })

# check:
# nrow(weightX)  # 10,354
# plot(colSums(weightX))
# hist(weightX[,1], breaks=100)
# abline(v=1, col="#AA0000")
# summary(weightX[,1])  #    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#                       # 0.0000  0.1480  0.4715  1.0000  1.1960 19.8400


# save(weightX, file = "weightX.Rdata")
# load("weightX.Rdata")


# ======  HELPER FUNCTIONS AND DEFINITIONS ===========================

# To weight the Jaccard Distance, we replace the sum and union
# functions with a sum over the expression weights



dWeightJaccard <- function(geneSetA, geneSetB, expWeights) {
  # A and B are gene sets, C is an expression set for
  # a cell line
  # We weight Jaccard Distance by the differential expression values
  num <- sum(expWeights[intersect(geneSetA, geneSetB)])
  denom <- sum(expWeights[union(geneSetA, geneSetB)])
  return(1 - (num / denom))
}


# We are repating code for FJD hash below - refactoring should
# merge this with the other function
WFJDhash <- hash()
WFJD <- function(A, B, C, drugs, pathN, wX, scal = FSCALE, reset=FALSE) {
  # This calculates an expression weighted Fuzzy Jaccard Distance
  # for drug target pathway neighborhoods.
  #
  # A: string: compound A
  # B: string: compound B
  # C: string: cell-line C
  #
  # drugs: a list of drug targets
  # pathN: a list of pathways that contain drug targets and their
  #        neighborhoods.
  # exp: weighted differential expressio values.
  #
  # This function maintains its results in a hash.
  # Call the function with WFJD(reset=TRUE) to reset the hash.
  #
  #
  if (reset) {
    WFJDhash <<- hash()
    return()
  }
  key <- sprintf("%s.%s.%s", A, B, C)
  dWFJ <- WFJDhash[[key]]
  if (is.null(dWFJ)) {  # not present in hash ...
    a <- pathN[[A]]
    b <- pathN[[B]]
    cellX <- wX[ , C]
    if (length(a$N0) == 0 | length(b$N0) == 0) {
      dWFJ <- 1
    } else {
      # calculate dFJ
      dWFJ <-        scal[1] * dWeightJaccard(a$N0, b$N0, cellX)
      dWFJ <- dWFJ + scal[2] * dWeightJaccard(a$N0, b$N1, cellX)
      dWFJ <- dWFJ + scal[3] * dWeightJaccard(a$N1, b$N0, cellX)
      dWFJ <- dWFJ + scal[4] * dWeightJaccard(a$N0, b$N2, cellX)
      dWFJ <- dWFJ + scal[5] * dWeightJaccard(a$N2, b$N0, cellX)
      dWFJ <- dWFJ + scal[6] * dWeightJaccard(a$N1, b$N1, cellX)
      dWFJ <- dWFJ + scal[7] * dWeightJaccard(a$N1, b$N2, cellX)
      dWFJ <- dWFJ + scal[8] * dWeightJaccard(a$N2, b$N1, cellX)
      dWFJ <- dWFJ + scal[9] * dWeightJaccard(a$N2, b$N2, cellX)
    }
    WFJDhash[[key]] <<- dWFJ      # store in hash
  }
  return(dWFJ)
}

A <- names(drugs)[1]
B <- names(drugs)[2]
C <- colnames(weightX)[1]


# Plot the Weighted Fuzzy Jaccard Distance for all drugs for one cell-line
# (Analytics only)
ptm <- proc.time() # Start the stopwatch...
N <- length(drugs)
dWFMat2 <- matrix(numeric(N * N), nrow=N)
cat("\n")
for (i in 1:N) {
  cat("=")
  for (j in i:N) {
    dWFMat2[i, j] <- dWFMat2[j, i] <- WFJD(names(drugs)[i],
                                           names(drugs)[j],
                                           colnames(weightX)[2],
                                           drugs, drugN, weightX)
  }
  if (! (i %% 20)) { cat("\n")}
}
cat("\n")
image(dWFMat2, col = colorRampPalette(c("#000000", "#5588FF", "#FFFFEE"))(12))
proc.time() - ptm  # How long did we take?
                   #    user  system elapsed
                   # 165.114  16.632 181.680
# That would be 3.8 hours  for all - unfortunately can't make that in time for the
# submission  :-(
# image(dWFMat1 - dWFMat2, col = colorRampPalette(c("#FF0000", "#000000", "#00FF00"))(12))
# plot(dMat, dFMat, xlim=c(0,1), ylim=c(0,1))



# [END]
