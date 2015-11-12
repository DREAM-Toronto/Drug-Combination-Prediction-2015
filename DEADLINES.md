# Challenge deadlines and format
There can be 3 submissions per round except for the Final Round.

## Subchallenge 1
There are three rounds plus the Final Round.

The submission format for 1A and 1B is the same: one zip file containing two csv files named prediction.csv, and combination_priority.csv.

For prediction.svn: three columns are CELL_LINE (normalized cell line name), COMBINATION_ID (two drug names separated by a dot) and PREDICTION (prediction value).

For combination_priority.csv: two columns are COMBINATION_ID (two drug names separated by a dot, the same as above) and CONFIDENCE (a decimal representation ranging from 0 to 1, the higher the more confident in the prediction regarding this combination).

### Passed
* Sept. 3, 2015: Challenges opened
* Sept. 10, 2015: Round 1 started
* Oct. 15, 2015: Submission started
* Nov. 6, 2015: Round 1 ended; Round 2 started

### Due
* **Dec. 17, 2015: Round 2 ends; Round 3 starts (updated on Nov 12))**
* Jan. 26, 2016: Round 3 ends; Final Round starts (1 submission only) (updated on Nov 12)
* Feb. 9, 2016: Final Round ends; Scoring starts (updated on Nov 12)
* Feb. 16, 2016: Scoring ends (updated on Nov 12)

## Subchallenge 2
There are two rounds plus the Final Round.

**Format**: One zipped File containing 'synergy_matrix.csv' and 'confidence_matrix.csv'

* synergy_matrix.csv
  * contain 85 cell lines(columns) by 370 drug combinations(rows)  ==>  this file is for leaderboard set
    * Prediction in this file should be either synergy(=1) or non-synergy(=0)
    * Non-synergy could be antagonism, no effect or additivity
  * Null predictions are not permitted
  * Column header (=first row in the csv file) must conain all 85 CELL_LINEs (note: first first column should be empty to form a matrix, e.g. ", cell_1, cell_2, ..., cell_85")
  * Row header (first column in the csv file) must contain all COMBINATION_IDs from the leader board, or final test set of challenge 2.
* confidence_matrix.csv
  * contain 85 cell lines(columns) by 370 drug combinations(rows)  ==>  this file is for final validation set
    * Continuous scores from 0 (no confidence) to 1 (most confidence) for each prediction
  * CONFIDENCE value cannot be null or NA.
  * Only valid CELL_LINE's, meaning controlled vocabulary cell line names (i.e. cell identifiers).
  * Only valid COMBINATION_ID's, meaning controlled vocabulary drug names separated by dot and alphabetically ordered.
  * List must contain 370 drug combinations for the leaderboard and 370 drug combinations for final submission
  * List must contain all 85 CELL_LINEs.

**Submission**<br />
1. Upload your prediction file to Synapse (under "Files" of an existing Synapse Project).<br />
2. Select the file in Synapse by opening the file's page.<br />
3. At the right hand side of the page click "Tools" > "Submit to Challenge".<br />
4. Select the Sub-challenge of interest.<br />
5. Provide a name for the team which created the results and a descriptive name for the submission.<br />
6. The confirmation dialog gives the address of the Synapse page where the leader board is displayed. Within 10 minutes your submission should appear in the corresponding leaderboard. If your submission is invalid, you will be notified via your registered e-mail.<br />

### Passed
* Sept. 3, 2015: Challenges opened
* Sept. 10, 2015: Round 1 started
* Oct. 15, 2015: Submission started

### Due
* **Dec 3, 2015: Round 1 ends; Round 2 starts (updated on Nov 12)**
* Jan. 26, 2016: Round 2 ends; Final Round starts (1 submission only) (updated on Nov 12)
* Feb. 9, 2016: Final Round ends; Scoring starts (updated on Nov 12)
* Feb. 16, 2016: Scoring ends (updated on Nov 12)
