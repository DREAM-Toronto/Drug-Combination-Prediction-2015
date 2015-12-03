# Retrieve PubChem data for all 120 drugs
# Should take around 2 mins 
# get drug compound names and smiles or cid 
raw <- read.csv(paste("../Challenge Data/Drug Synergy Data/", "Drug_info_release.csv", sep=""),
                head=FALSE,
                stringsAsFactors=FALSE)

DESTFOLDER = "../ExternalData/PubChemData/"

# for each row, query and save PubChem data using PUG REST
for (i in 2:nrow(raw)) {

  compound = raw[i,1]
  queries = raw[i,7]
  
  if(queries != "") {
    
    # some 'drugs' consist of more than one compounds
    # with SMILES or CID separated by ";"
    splitQueries = strsplit(queries, ";")[[1]]
    
    for (i in 1:length(splitQueries)) {
      
      query = splitQueries[i]
      
      # if query contains only digits, it is a CID, else it is a SMILES
      if (regexpr("[0-9]+", query) == TRUE) {
        
        # query is a CID
        inputQuery = paste("/cid/", query, sep = '')
        
      } else {
        
        # query is a SMILES 
        # SMILES contain special characters that need to be URL-encoded
        # For example:
        #     SMILES = c1cc(ccc1[C@H](CCO)NC(=O)C2(CCN(CC2)c3c4cc[nH]c4ncn3)N)Cl
        #     URL ENCODED SMILES = c1cc%28ccc1%5BC%40H%5D%28CCO%29NC%28%3DO%29C2%28CCN%28CC2%29c3c4cc%5BnH%5Dc4ncn3%29N%29Cl
        smiles = URLencode(query, reserved = TRUE)    
        inputQuery = paste("/smiles/", smiles, sep = '')
        
      }
      
      # Prepare the URL for HTTP
      prolog = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
      inputDomain = "/compound"
      output = "/JSON"
      operations = "/?response_type=display" # Operation Options "/?response_type=display" 
      input = paste(inputDomain, inputQuery, sep = '')
      
      url = paste(prolog, input, output, operations, sep = '')
    
      # compounds of the same drug will be distinguisehd by (compound #) 
      if (length(splitQueries) != 1) {
        fname = paste(compound, "(compound", i, ")", sep = '')
      } else {
        fname = compound
      }
      
      # destination file to save the json
      file = paste(DESTFOLDER, fname, ".json", sep = '')
      
      # download json output into the destination file 
      download.file(url, destfile = file, mode = "w")
    }
  } 
  
}