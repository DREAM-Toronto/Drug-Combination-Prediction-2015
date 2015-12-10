# Requires db and user be set up first
# Database name: dream2015
# Recommend downloading MySQL workbench https://dev.mysql.com/downloads/workbench/
DB <- "dream2015"
USER <- "ehsueh"
PW <- "dream2015"
HOST <- "localhost"
DATADIR <- "../Challenge Data/Drug Synergy Data/"

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
query_create_tables <- "create table if not exists drug_synergy_data (
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
);"

dbSendQuery(connection, query_create_tables)


# ==========================================================
# LOADING DATA
# ==========================================================

# load raw data
fileList <- c("ch1_leaderBoard_monoTherapy.csv",
"ch1_test_monoTherapy.csv",
"ch1_train_combination_and_monoTherapy.csv",
"ch2_leaderBoard_monoTherapy.csv",
"ch2_test_monoTherapy.csv")

for (file in fileList) {
  
  raw <- read.csv(paste(DATADIR, file, sep=""),
                  head = TRUE,
                  stringsAsFactors = FALSE,
                  sep = ",")

  # source column
  src <- rep(0,nrow(raw))
  src[1:nrow(raw)] <- file
  
  # convert into a data frame
  df <- data.frame(src, raw);
  
  # load into db
  dbWriteTable(connection, value = df, name = "drug_synergy_data", append = TRUE, row.names = FALSE)
  
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

# extract IC50, H, Einf features and load into feature tables
drug1_features <- getQueryResults("select compound_a, cell_line, ic50_a, h_a, einf_a from drug_synergy_data")
drug2_features <- getQueryResults("select compound_b, cell_line, ic50_b, h_b, einf_b from drug_synergy_data")

dbDisconnect(connection)
# ==========================================================
# EXIT and CLOSE CONNECTION
# ==========================================================
on.exit(dbDisconnect(connection))