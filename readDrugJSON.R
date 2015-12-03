setwd("~/Desktop/JSON_Practice")

if (! require(jsonlite, quietly=TRUE)) {
  # for plotting dose-response surfaces
  install.packages("jsonlite",
                   repos="http://cran.us.r-project.org")
  library(jsonlite)
}

readDrugJSON <- function(fn){
  if (! require(jsonlite, quietly=TRUE)) {
    # for plotting dose-response surfaces
    install.packages("jsonlite",
                     repos="http://cran.us.r-project.org")
    library(jsonlite)
  }
  
  x <- fromJSON(fn,FALSE)
  chemPhyProperties <- x$Record$Section[[4]]
  computedProperties <- chemPhyProperties$Section[[1]]
  # experimentalProperties <- chemPhyProperties$Section[[2]]
  
  computedPropertyObj <- list();
  
  for(i in 1:length(computedProperties$Section)){
    curProperty <- computedProperties$Section[[i]]$Information[[1]]
    newPropertyObj <- {}
    curPropertyName <-curProperty[2][[1]]
    # newPropertyObj$Name <- curPropertyName
    newPropertyObj$Value <- curProperty[3][[1]]
    if(length(curProperty)>3){
      newPropertyObj$Unit <- curProperty[4][[1]]
    }
    # computedPropertyObj[[i]] <- newPropertyObj
    computedPropertyObj[[curPropertyName]] <- newPropertyObj
  }
  
  return(computedPropertyObj)
}

a <- readDrugJSON("CID_9444.json")

# x <- fromJSON("CID_9444.json",FALSE)
# 
# chemPhyProperties <- x$Record$Section[[4]]
# computedProperties <- chemPhyProperties$Section[[1]]
# # experimentalProperties <- chemPhyProperties$Section[[2]]
# 
# computedPropertyObj <- list();
# 
# for(i in 1:length(computedProperties$Section)){
#   curProperty <- computedProperties$Section[[i]]$Information[[1]]
#   newPropertyObj <- {}
#   curPropertyName <-curProperty[2][[1]]
#   # newPropertyObj$Name <- curPropertyName
#   newPropertyObj$Value <- curProperty[3][[1]]
#   if(length(curProperty)>3){
#     newPropertyObj$Unit <- curProperty[4][[1]]
#   }
#   # computedPropertyObj[[i]] <- newPropertyObj
#   computedPropertyObj[[curPropertyName]] <- newPropertyObj
# }

