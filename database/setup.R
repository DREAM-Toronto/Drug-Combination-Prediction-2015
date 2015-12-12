# Requires db and user be set up first
# Database name: dream2015
# Recommend downloading MySQL workbench https://dev.mysql.com/downloads/workbench/
DREAMDIR <- "/media/ehsueh/Data/projects/dream/src/Drug-Combination-Prediction-2015/" # if you didn't source your RStudio
setwd(DREAMDIR)

DB <- "dream2015"
USER <- "ehsueh"
PW <- "dream2015"
HOST <- "localhost"
DATADIR <- paste("../Challenge Data/", sep = "")

# ==========================================================
# RESOLVE DEPENDENCIES
# ==========================================================
# might need to execute this in the terminal first
# sudo apt-get install r-cran-rmysql
if (! require(RMySQL, quietly=TRUE)) {
  # for plotting dose-response surfaces
  install.packages("RMySQL",
                   repos="http://cran.us.r-project.org")
  library(RMySQL)
}

# ==========================================================
# CREATING A DATABASE BY RUNNING setup.sql
# ==========================================================
# system(paste("mysql -u ", USER, " -D ", DB, " -p ", PW, " < ./database/setup.sql", sep = ""))

# ==========================================================
# OPENING A CONNECTION
# ==========================================================
# opening db connection
connection <- dbConnect(MySQL(),
                        user = USER,
                        password = PW,
                        dbname = DB,
                        host = HOST)

# ==========================================================
# CREATING TABLES
# ==========================================================
# contains drug synergy data from provided ch*.csv files  

query_create_table <- "

create table if not exists raw_drug_info (
  id int not null auto_increment primary key,
  challenge_name varchar(25) not null,
  target varchar(50) not null,
  hba int,
  clogp decimal(5,3),
  hbd int,
  lipinski int,
  smiles varchar(300) not null,
  pubchem_id int,
  mw decimal(5,3)
);
"
dbSendQuery(connection, query_create_table)

query_create_table <- "

create table if not exists raw_cell_info (
  id int not null auto_increment primary key,
  sanger_name varchar(25) not null,
  ccle_name varchar(50) not null,
  alt_name varchar(50) not null,
  disease_area varchar(50) not null,
  tissue_general varchar(50) not null,
  cosmic int not null
);

"
dbSendQuery(connection, query_create_table)

query_create_table <- "

create table if not exists raw_drug_synergy_data (
  id int not null auto_increment primary key,
  src varchar(30) not null,
  cell_line varchar(30) not null,
  compound_a varchar(30) not null,
  compound_b varchar(30) not null,
  max_conc_a decimal(5,2) not null, # precision 5, scale 2 (i.e. xxx.xx)
  max_conc_b decimal(5,2) not null,
  ic50_a decimal(20,10) not null,
  ic50_b decimal(20,10) not null,
  h_a decimal(20,10) not null,
  h_b decimal(20,10) not null,
  einf_a decimal(20,10) not null,
  einf_b decimal(20,10) not null,
  synergy_score decimal(9,6) not null,
  qa int not null,
  combination_id varchar(30) not null
);

"
dbSendQuery(connection, query_create_table)

# ==========================================================
# LOADING DATA
# ==========================================================

# load raw cell and drug info csv
z <- c("drug", "cell")
fileList <- c("Drug Synergy Data/Drug_info_release.csv", "Sanger Molecular Data/cell_info.csv")

processComponentName <- function(compo, file) {
  tableName1 <- paste("raw", compo, "info", sep = "_")
  tableName2 <- paste(compo, "name", sep = "_")
  raw <- read.csv(paste(DATADIR, file , sep=""),
                  head = TRUE,
                  stringsAsFactors = FALSE,
                  sep = ",")
  df <- data.frame(raw);
  dbWriteTable(connection, value = df, name = tableName1, append = TRUE, row.names = FALSE)
  df2 <- df[1]
  names(df2)[1] <- "name"
  dbWriteTable(connection, value = df2, name = tableName2, append = TRUE, row.names = FALSE)
  return(df2)
}

dfDrug <- processComponentName(z[1], fileList[1])
dfCell <- processComponentName(z[2], fileList[2])

# load raw monotherapy csv
fileList <- c("ch1_leaderBoard_monoTherapy.csv",
"ch1_test_monoTherapy.csv",
"ch1_train_combination_and_monoTherapy.csv",
"ch2_leaderBoard_monoTherapy.csv",
"ch2_test_monoTherapy.csv")

for (file in fileList) {
  
  raw <- read.csv(paste(DATADIR, "Drug Synergy Data/", file, sep=""),
                  head = TRUE,
                  stringsAsFactors = FALSE,
                  sep = ",")

  # source column
  src <- rep(0,nrow(raw))
  src[1:nrow(raw)] <- file
  
  # convert into a data frame
  df <- data.frame(src, raw);
  
  # load into db
  dbWriteTable(connection, value = df, name = "raw_drug_synergy_data", append = TRUE, row.names = FALSE)
  
}

# helper function for getting query results
getQueryResults <- function(query) {
  resultSet <- dbSendQuery(connection, query)
  results <- fetch(res = resultSet, n = -1) # -1 means get all rows
  if (dbHasCompleted(resultSet)) { # must close the result set before continuing
    dbClearResult(resultSet)
  }
  return (results)
}

# type can be "drug, cell, drug_drug, drug_cell, or drug_drug_cell"
addFeature <- function(feature, type) {
  if (type %in% c("drug", "cell", "drug_drug", "drug_cell", "drug_drug_cell")) {
    # insert if unique
    # else return the id of specified feature
    # TO-DO's: Not a priority, but it will be nice if the id does not increment when insertions failed.
    #          Right now, inserting duplicate values results in id incrementation. 
    #          I.e. we are left with ugly id's that jumps around.
    query <- paste("insert into ", type, "_feature(feature) values('", feature,
                   "') on duplicate key update id = last_insert_id(id), feature ='", feature, "';", sep = "")
    query1 <- "select last_insert_id();"
    dbSendQuery(connection, query)
    resultSet <- dbSendQuery(connection, query1)
    id <- fetch(resultSet, n = 1)
    dbClearResult(resultSet)
    return(id)
  } else {
    # maybe we should crash here
    print("Incorrect feature type.")
  }
}

# extract IC50, H, Einf features and load into feature tables
drug1_features <- getQueryResults("select d.id as drug_id, c.id as cell_id,
                                   avg(r.ic50_a/r.max_conc_a) as ic50,
                                   avg(r.h_a*r.max_conc_a) as h,
                                   avg(r.einf_a) as einf
                                   from raw_drug_synergy_data r, drug_name d, cell_name c 
                                   where r.compound_a = d.name and r.cell_line = c.name
                                   group by r.compound_a, r.cell_line
                                   order by r.compound_a, r.cell_line asc;")
drug2_features <- getQueryResults("select d.id as drug_id, c.id as cell_id,
                                   avg(r.ic50_b/r.max_conc_b) as ic50,
                                   avg(r.h_b*r.max_conc_b) as h,
                                   avg(r.einf_b) as einf
                                   from raw_drug_synergy_data r, drug_name d, cell_name c 
                                   where r.compound_b = d.name and r.cell_line = c.name
                                   group by r.compound_b, r.cell_line
                                   order by r.compound_b, r.cell_line asc;")

drug_features <- rbind(drug1_features, drug2_features)

params <- c("ic50", "h", "einf")
# insert ic50_avg_normalized, h_avg_normalized, einf_avg_normalized <-- only one set per (drug, cell) pair

for (i in seq(1:3)) {
  id <- addFeature(paste(params[i], "_avg_normalized", sep = ""), "drug_cell")
  df_avg_n <- data.frame(drug_features[,1:2], id, drug_features[,2+i])
  names(df_avg_n)[3] <- "feature_id"
  names(df_avg_n)[4] <- "value"
  dbWriteTable(connection, value = df_avg_n, name = "drug_cell", append = TRUE, row.names = FALSE)
}

# TO-DO's: 
# insert different versions of hill function parameters for example:
#   insert ic50_raw, h_raw, einf_raw <-- there could be multiple sets of values per (drug, cell) pair
#   insert ic50_median, h_median, einf_median <-- only one set per (drug, cell) pair
#   insert ic50_median_normalized, h_median_normalized, einf_median_normalized <-- only one set per (drug, cell) pair

# ==========================================================
# EXIT and CLOSE CONNECTION
# ==========================================================
on.exit(dbDisconnect(connection))