# DREAM Challenge Toronto 2015
## Drug Combination Prediction

### Data

#### From DREAM

Dose response curves
* per compound
* per cell line
* for combinations (provide training set for learning algorithm)

Molecules
* molecular weight
* solubility
* chemical fragments
* shape
* QSAR

Genomic
* CNVs (copy number variations)
    * parts of genome are {duplicated|lost}
    * changes in expression levels -> gene dosage effects
    * most important for gene loss
* SNPs - structural, splice sites, regulation

Methylation
* histones -> expression

#### External Data

Epistatic interactions
* A or B by itself not effective, need to stop both
* Pathway C might regulate/interact with A/B

Network data
* pathways (priority to pathways involved in cancer?)
    * cell growth
    * genome stability
    * loss of contact inhibition
    * metastasis
* PPi
* co-expression data (quantitative relationships b/w proteins and genes)

### Sub-Challenge 1

Feature -> ML -> Combination scores

### Sub-Challenge 2

**Round 1 DeadLine has been changed to Dec 3.**(See DeadLine File for more details)<br />
1. Infer Drug Synergy without experimental synergy score.<br />
2. Making drug predictions based on prior knowledge.<br />
3. We are allowed to use ** molecular data for the cell lines, cell response data for all respective    mono-therapies, chemical information and putative targets of the compounds**<br />
4. Require 740 drug combinations without overlapping with Subchallenge 1. <br />
   -- 4a. A leaderboard set (370 combinations) and a final validation set (370 combinations) are required.<br />
   -- 4b. A full synergy-prediction-matrix (Score either 1 or 0)<br />
   -- 4c. A full synergy-confidence-matrix (Score ranging from 0 to 1)<br />
5. The above data will be used to score accuracy of predictions.<br />
6. Justify our stratifications and synergy predictions, by giving an algorithm and rationally translatable as biomarkers.<br />


## Checking in...

* Boris (Hyginn)
* Emma (ehsueh)
* Julian (thejmazz)
* Nathan (njia95)
* Pruthvi (Pruthv1)
* Ashley (AshleyWWW)
* Jack (c5chenpe)
* Bhawan (B-1P)

## Groups & Responsibilities

### Workflow

* Targets: Nathan, Jack
* Data Preprocessing: Ashley, Fred
* Data Reduction: Boris
* Cross Validation: Chris, Julian
* Drug Similarity: Bhawan, Zach, Ricardo, Emma
* Cell Line Similarity: Pruthvi, Jenny

### Support
* Scoring:
* Submission: Jack, Nathan
* Documentation: Julian
* QA:
